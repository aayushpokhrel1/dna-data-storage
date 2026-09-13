"""M3 recovery benchmark: encode -> corrupt -> decode, measured vs error rate.

Builds on the M2 ECC layer (`ecc.encode`/`ecc.decode`) and the M3 channel
(`channel.corrupt`). The deliverable is the recovery-vs-rate sweep; error rates
are parameters, not baked-in constants. No rate here is a literature value
(PENDING CITATION, see the M3 spec).

Metrics (definitions fixed in-house):
  - recovery at a rate = fraction of trials that decode exactly to the input;
  - density = payload bits / total encoded nt (bits/nt);
  - code rate = payload bytes / encoded record bytes (parity + CRC + index all
    count as overhead; the formula is documented at `code_rate`).
"""
import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
import random

import channel
import codec
import constraints
import ecc
import fountain
import markercode
import screen_codec
from reedsolo import ReedSolomonError

RESULTS_DIR = os.path.join(os.path.dirname(__file__), os.pardir, "results")
RESULTS_PATH = os.path.join(RESULTS_DIR, "benchmark.json")

# Cited per-base error context (verified 2026-09-11) from the review "Uncertainties
# in synthetic DNA-based data storage" (PMC8191772): phosphoramidite synthesis ~0.7%
# (substitution 0.5%, insertion 0.1%, deletion 0.1%); Illumina sequencing <1%;
# nanopore ~10%. Used to set the sub:ins:del ratio and to frame the swept range; the
# numbers are cited in paper/references.bib. See docs/M3 spec.
LIT_RATIO = (5, 1, 1)  # substitution:insertion:deletion ~ 0.5:0.1:0.1 (PMC8191772)


def code_rate(data_len, oligos, meta):
    """Density and overhead of an encoded pool.

    bits_per_nt = data_len * 8 / total_nt, total_nt = sum of oligo lengths.

    code_rate = data_len / (data_len + parity_bytes_equivalent) reduces to
    data_len / (n_records * data_bytes): the outer RS block rate, counting the
    parity records and the zero-padding of the last block. It does NOT count the
    per-record CRC, the oligo index, or the base-3 rotating expansion; those are
    all captured end-to-end by bits_per_nt, which is the true density metric.
    n_records = sum(block_sizes) + parity_records * n_blocks.
    """
    e = meta["ecc"]
    D, P = e["data_bytes"], e["parity_records"]
    total_nt = sum(len(o) for o in oligos)
    n_blocks = len(e["block_sizes"])
    n_records = sum(e["block_sizes"]) + P * n_blocks
    parity_bytes_equivalent = n_records * D - data_len
    denom = data_len + parity_bytes_equivalent
    return {
        "bits_per_nt": data_len * 8 / total_nt if total_nt else 0.0,
        "code_rate": data_len / denom if denom else 0.0,
        "n_oligos": len(oligos),
        "total_nt": total_nt,
    }


def constraint_stats(oligos):
    """Per-oligo GC (min/mean/max) and per-oligo max homopolymer run.

    Each oligo is one molecule, so runs are measured per oligo and the reported
    run is the max over oligos. Measuring over joined oligos would invent runs
    at the boundaries.
    """
    gcs = [constraints.gc_fraction(o) for o in oligos]
    runs = [constraints.max_homopolymer_run(o) for o in oligos]
    return {
        "gc_min": min(gcs, default=0.0),
        "gc_mean": sum(gcs) / len(gcs) if gcs else 0.0,
        "gc_max": max(gcs, default=0.0),
        "max_homopolymer_run": max(runs, default=0),
    }


def recovery_at_rate(data, rate, trials, seed, ratio=(1, 1, 1), coverage=1, **ecc_kwargs):
    """Fraction of `trials` encode->corrupt->decode passes that recover `data`.

    `rate` is the per-base error rate, split across (sub, ins, del) by `ratio`.
    `coverage` is the number of sequencing reads per oligo. Each trial uses a
    distinct channel seed derived from `seed`. A RecoveryError or a mismatch
    counts as a failure and never propagates.
    """
    total = sum(ratio)
    p_sub, p_ins, p_del = (rate * r / total for r in ratio)
    oligos, meta = ecc.encode(data, **ecc_kwargs)
    ok = 0
    for t in range(trials):
        reads = channel.corrupt(oligos, p_sub=p_sub, p_ins=p_ins, p_del=p_del,
                                coverage=coverage, seed=seed + t)
        try:
            if ecc.decode(reads, meta) == data:
                ok += 1
        except (ecc.RecoveryError, ReedSolomonError):
            pass  # both are decode-failure modes -> counts as non-recovery
    return ok / trials if trials else 0.0


def density_sweep(payload_sizes, n_data_records=32, parity_records=4):
    """Density vs `payload_bytes`, with the data-record count held fixed.

    Fixing `n_data_records` (and putting them in one block) keeps the parity
    fraction constant, so the trend isolates how per-oligo overhead (the index and
    the CRC) amortizes as the payload grows: bits/nt rises toward the rotating
    ceiling (log2 3 ~ 1.585) scaled by the code rate. Payload is deterministic.
    """
    out = {"payload_bytes": [], "bits_per_nt": [], "code_rate": []}
    for pb in payload_sizes:
        data = bytes((i * 7 + 1) % 256 for i in range(pb * n_data_records))
        oligos, meta = ecc.encode(data, data_bytes=pb, parity_records=parity_records,
                                  block_records=n_data_records)
        cr = code_rate(len(data), oligos, meta)
        out["payload_bytes"].append(pb)
        out["bits_per_nt"].append(cr["bits_per_nt"])
        out["code_rate"].append(cr["code_rate"])
    return out


def sweep(data, rates, trials, seed, **kwargs):
    """Run `recovery_at_rate` over `rates`; return aligned rate/recovery lists."""
    recovery = [recovery_at_rate(data, r, trials, seed, **kwargs) for r in rates]
    config = {"trials": trials, "seed": seed, **kwargs}
    return {"rates": list(rates), "recovery": recovery, "config": config}


def _family_recovery(data, rate, trials, seed, base, layer, enc_kwargs,
                     coverage, ratio=LIT_RATIO, p_drop=0.0):
    """Fraction of trials that a given (base codec, ECC layer) family recovers."""
    total = sum(ratio)
    p_sub, p_ins, p_del = (rate * r / total for r in ratio)
    oligos, meta = layer.encode(data, base=base, **enc_kwargs)
    ok = 0
    for t in range(trials):
        reads = channel.corrupt(oligos, p_sub=p_sub, p_ins=p_ins, p_del=p_del,
                                p_drop=p_drop, coverage=coverage, seed=seed + t)
        try:
            if layer.decode(reads, meta, base=base) == data:
                ok += 1
        except (ecc.RecoveryError, fountain.RecoveryError, ReedSolomonError):
            pass
    return ok / trials if trials else 0.0


def compare_families(data, rates, overhead=1.0, coverage=6, trials=8,
                     data_bytes=16, seed=0):
    """Density and recovery for all four {base codec} x {ECC layer} families,
    at a matched redundancy (`overhead`) and coverage. This is the cross-family
    comparison the paper's Pareto frontier is built from.

    `overhead` is the extra fraction of records: RS parity_records = overhead * K,
    fountain droplets = (1 + overhead) * K, so both families carry the same
    redundancy and only the codec/ECC family differs.
    """
    K = max(1, -(-len(data) // data_bytes))
    parity = max(1, round(K * overhead))
    rs = dict(data_bytes=data_bytes, parity_records=parity, block_records=K)
    ft = dict(data_bytes=data_bytes, overhead=overhead)
    families = {
        "rotating+RS": (codec, ecc, rs),
        "rotating+fountain": (codec, fountain, ft),
        "screening+RS": (screen_codec, ecc, rs),
        "screening+fountain": (screen_codec, fountain, ft),
    }
    out = {}
    for name, (base, layer, enc) in families.items():
        oligos, _ = layer.encode(data, base=base, **enc)
        bits_per_nt = len(data) * 8 / sum(len(o) for o in oligos)
        recovery = [_family_recovery(data, r, trials, seed, base, layer, enc, coverage)
                    for r in rates]
        threshold = max([r for r, v in zip(rates, recovery) if v >= 1.0], default=0.0)
        out[name] = {"bits_per_nt": bits_per_nt, "recovery": recovery,
                     "threshold": threshold}
    return {"rates": list(rates), "overhead": overhead, "coverage": coverage,
            "families": out}


def _enc_kwargs(layer, ov, K, data_bytes, block_records):
    """Build encode kwargs for a family at redundancy `ov` (extra fraction).

    RS parity is per block, so it scales with the block size (not global K); that
    keeps the redundancy fraction ov matched to fountain and within the GF(256)
    block limit when K forces multiple blocks.
    """
    if layer is ecc:
        bsize = block_records or K
        return dict(data_bytes=data_bytes, parity_records=max(1, round(bsize * ov)),
                    block_records=bsize)
    return dict(data_bytes=data_bytes, overhead=ov)


def cost_grid(data, rate, coverages, overheads, base=codec, layer=ecc,
              data_bytes=16, trials=8, seed=0, p_drop=0.0, block_records=None):
    """Recovery over a (coverage x redundancy) grid at a fixed error rate.

    Sequencing coverage costs reads; redundancy costs synthesized bases. The grid,
    and the iso-recovery frontier through it, show the cheapest budget that still
    recovers. Rows are coverages, columns are overheads.
    """
    K = max(1, -(-len(data) // data_bytes))
    recovery = []
    for cov in coverages:
        row = []
        for ov in overheads:
            enc = _enc_kwargs(layer, ov, K, data_bytes, block_records)
            row.append(_family_recovery(data, rate, trials, seed, base, layer, enc,
                                        cov, p_drop=p_drop))
        recovery.append(row)
    return {"coverages": list(coverages), "overheads": list(overheads),
            "rate": rate, "p_drop": p_drop, "recovery": recovery}


def dropout_study(data, overheads, p_drop=0.3, coverage=4, data_bytes=16,
                  block_records=64, base=codec, trials=8, seed=0):
    """RS vs fountain recovery vs redundancy at high whole-oligo dropout.

    Uses a large K forced into multiple RS blocks (`block_records` < K), the regime
    where RS's per-block erasure budget is fragile to uneven dropout while a fountain
    code pools droplets globally. Reports whatever it shows, honestly.
    """
    K = max(1, -(-len(data) // data_bytes))
    out = {"overheads": list(overheads), "p_drop": p_drop, "coverage": coverage,
           "k": K, "block_records": block_records}
    for name, layer in (("RS", ecc), ("fountain", fountain)):
        recs = []
        for ov in overheads:
            enc = _enc_kwargs(layer, ov, K, data_bytes, block_records)
            recs.append(_family_recovery(data, 0.0, trials, seed, base, layer, enc,
                                         coverage, p_drop=p_drop))
        out[name] = recs
    return out


def _rotating_bases(n, seed):
    rng = random.Random(seed)
    out, prev = [], "A"
    for _ in range(n):
        b = rng.choice([x for x in "ACGT" if x != prev])
        out.append(b)
        prev = b
    return "".join(out)


def _correct_fraction(a, b):
    n = min(len(a), len(b))
    return sum(x == y for x, y in zip(a[:n], b[:n])) / n if n else 0.0


def indel_confinement(bases_len=512, del_rates=(0.0, 0.005, 0.01, 0.02, 0.05),
                      period=16, trials=10, seed=0):
    """Fraction of bases recovered vs deletion rate, with marker resync vs without.

    Shows the S6 result: a deletion without markers shifts every downstream base
    (recovery falls off fast), while the marker layer re-anchors and confines the
    damage to one run (recovery stays high).
    """
    with_markers, without = [], []
    for r in del_rates:
        wm, wo = [], []
        for t in range(trials):
            bases = _rotating_bases(bases_len, seed + t)
            (read_m,) = channel.corrupt([markercode.encode(bases, period=period)],
                                        p_del=r, seed=seed + t)
            wm.append(_correct_fraction(
                markercode.decode(read_m, bases_len, period=period), bases))
            (read_o,) = channel.corrupt([bases], p_del=r, seed=seed + t)
            wo.append(_correct_fraction(read_o, bases))
        with_markers.append(sum(wm) / len(wm))
        without.append(sum(wo) / len(wo))
    return {"del_rates": list(del_rates), "with_markers": with_markers,
            "without": without, "period": period}


def default_family_comparison():
    """The standard cross-family run, single source for the JSON and the figure."""
    data = bytes(range(256)) * 2  # 512 bytes
    rates = [0.0, 0.005, 0.01, 0.02, 0.05, 0.1]
    return compare_families(data, rates, overhead=1.0, coverage=8, trials=8,
                            data_bytes=16)


def default_cost_study():
    """The standard S4 run: the coverage x redundancy grid (for the frontier-winning
    screening+RS family) and the large-K dropout RS-vs-fountain study. Single source
    for the JSON and the figures."""
    grid = cost_grid(bytes(range(256)) * 2, rate=0.02, coverages=[1, 2, 4, 8, 16],
                     overheads=[0.25, 0.5, 1.0, 2.0], base=screen_codec,
                     trials=8, data_bytes=16)
    big = bytes((i * 131 + 7) % 256 for i in range(3000))  # K ~ 188, multi-block RS
    drop = dropout_study(big, overheads=[0.25, 0.5, 0.75, 1.0, 1.5], p_drop=0.3,
                         coverage=1, block_records=64, trials=6)
    return {"grid": grid, "dropout": drop}


def main(smoke=False):
    """Run a small sweep, write results/benchmark.json, print a summary."""
    if smoke:
        data = bytes(range(16))
        rates = [0.0, 0.01, 0.05]
        trials = 3
        coverage = 10
        ecc_kwargs = {"data_bytes": 8, "parity_records": 4}
    else:
        data = bytes(range(256))
        # spans synthesis-only (~0.7%) through nanopore-scale (~10%) per-base error
        rates = [0.0, 0.005, 0.01, 0.02, 0.05, 0.1]
        trials = 20
        coverage = 10
        ecc_kwargs = {"data_bytes": 8, "parity_records": 4}

    result = sweep(data, rates, trials, seed=0, ratio=LIT_RATIO,
                   coverage=coverage, **ecc_kwargs)
    oligos, meta = ecc.encode(data, **ecc_kwargs)
    payload = {
        "sweep": result,
        "code_rate": code_rate(len(data), oligos, meta),
        "constraint_stats": constraint_stats(oligos),
        "literature_context": {
            "sub_ins_del_ratio": LIT_RATIO,
            "coverage": coverage,
            "note": "per-base error range and sub:ins:del ratio framed on PMC8191772 "
                    "(synthesis ~0.7%: sub 0.5 / ins 0.1 / del 0.1; Illumina <1%; "
                    "nanopore ~10%); cited in paper/references.bib",
        },
    }
    os.makedirs(RESULTS_DIR, exist_ok=True)
    with open(RESULTS_PATH, "w") as f:
        json.dump(payload, f, indent=2)

    cr = payload["code_rate"]
    cs = payload["constraint_stats"]
    print(f"data    : {len(data)} bytes -> {cr['n_oligos']} oligos, "
          f"{cr['total_nt']} nt")
    print(f"density : {cr['bits_per_nt']:.3f} bits/nt   code rate: {cr['code_rate']:.3f}")
    print(f"GC      : {cs['gc_min']:.3f}/{cs['gc_mean']:.3f}/{cs['gc_max']:.3f} "
          f"(min/mean/max)   max run (per oligo): {cs['max_homopolymer_run']}")
    for r, rec in zip(result["rates"], result["recovery"]):
        print(f"  rate {r:.3f} -> recovery {rec:.2f}")
    print(f"wrote   : {os.path.normpath(RESULTS_PATH)}")

    if not smoke:  # cross-family comparison (S3)
        fam = default_family_comparison()
        fam_path = os.path.join(RESULTS_DIR, "family_comparison.json")
        with open(fam_path, "w") as f:
            json.dump(fam, f, indent=2)
        print("\ncross-family (overhead {:.1f}x, coverage {}):".format(
            fam["overhead"], fam["coverage"]))
        for name, d in fam["families"].items():
            print(f"  {name:20s} density={d['bits_per_nt']:.3f} bits/nt   "
                  f"full-recovery threshold={d['threshold']:.3f}")
        print(f"wrote   : {os.path.normpath(fam_path)}")

        cost = default_cost_study()  # S4 cost study
        cost_path = os.path.join(RESULTS_DIR, "cost_study.json")
        with open(cost_path, "w") as f:
            json.dump(cost, f, indent=2)
        dr = cost["dropout"]
        print(f"\ndropout study (K={dr['k']}, blocks {dr['block_records']}, "
              f"p_drop {dr['p_drop']}, overheads {dr['overheads']}):")
        print(f"  RS       {['%.2f' % v for v in dr['RS']]}")
        print(f"  fountain {['%.2f' % v for v in dr['fountain']]}")
        print(f"wrote   : {os.path.normpath(cost_path)}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="M3 recovery benchmark")
    parser.add_argument("--smoke", action="store_true",
                        help="run a tiny, fast sweep")
    args = parser.parse_args()
    main(smoke=args.smoke)
