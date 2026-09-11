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
import channel
import constraints
import ecc
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


def sweep(data, rates, trials, seed, **kwargs):
    """Run `recovery_at_rate` over `rates`; return aligned rate/recovery lists."""
    recovery = [recovery_at_rate(data, r, trials, seed, **kwargs) for r in rates]
    config = {"trials": trials, "seed": seed, **kwargs}
    return {"rates": list(rates), "recovery": recovery, "config": config}


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


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="M3 recovery benchmark")
    parser.add_argument("--smoke", action="store_true",
                        help="run a tiny, fast sweep")
    args = parser.parse_args()
    main(smoke=args.smoke)
