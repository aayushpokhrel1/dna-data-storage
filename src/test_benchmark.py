"""M3 self-checks for the recovery benchmark. Run: python src/test_benchmark.py"""
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
import benchmark
import ecc

DATA = bytes(range(64))


def test_zero_error_full_recovery():
    r = benchmark.recovery_at_rate(DATA, 0.0, trials=5, seed=0, parity_records=4)
    assert r == 1.0, f"expected full recovery at rate 0, got {r}"


def test_low_rate_within_budget_recovers():
    # generous parity, tiny rate -> should always recover (deterministic seed)
    r = benchmark.recovery_at_rate(DATA, 0.004, trials=5, seed=0,
                                   data_bytes=8, parity_records=8, block_records=64)
    assert r == 1.0, f"expected full recovery within budget, got {r}"


def test_coverage_improves_recovery():
    # At a per-base rate that wrecks most oligos with a single read, sequencing
    # coverage (multiple reads, keep a clean one) restores recovery.
    lo = benchmark.recovery_at_rate(DATA, 0.01, trials=3, seed=0,
                                    coverage=1, data_bytes=8, parity_records=4)
    hi = benchmark.recovery_at_rate(DATA, 0.01, trials=3, seed=0,
                                    coverage=25, data_bytes=8, parity_records=4)
    assert lo < 1.0, f"expected single-read recovery to fail at 1%, got {lo}"
    assert hi == 1.0, f"expected coverage to restore recovery, got {hi}"


def test_high_rate_degrades_without_crashing():
    r = benchmark.recovery_at_rate(DATA, 0.9, trials=5, seed=0, parity_records=2)
    assert 0.0 <= r < 1.0, f"expected degraded recovery at high rate, got {r}"


def test_sweep_shape():
    res = benchmark.sweep(DATA, [0.0, 0.5], trials=3, seed=0, parity_records=4)
    assert res["rates"] == [0.0, 0.5]
    assert len(res["recovery"]) == 2
    assert res["recovery"][0] == 1.0


def test_density_sweep_rises_toward_ceiling():
    import math
    sweep = benchmark.density_sweep([8, 16, 32, 64, 128], n_data_records=32,
                                    parity_records=4)
    bpn = sweep["bits_per_nt"]
    ceiling = math.log2(3)  # rotating code payload ceiling ~1.585 bits/nt
    assert bpn == sorted(bpn), f"density must rise with payload_bytes, got {bpn}"
    assert all(0 < b < ceiling for b in bpn), f"density out of (0, log2 3): {bpn}"
    assert bpn[-1] > bpn[0], "large payload must amortize index/CRC overhead"


def test_compare_families_reports_all_four():
    data = bytes(range(256))
    res = benchmark.compare_families(data, [0.0, 0.02], overhead=1.0, coverage=6, trials=3)
    fam = res["families"]
    assert set(fam) == {"rotating+RS", "rotating+fountain",
                        "screening+RS", "screening+fountain"}
    for name, d in fam.items():
        assert d["bits_per_nt"] > 0, name
        assert len(d["recovery"]) == 2, name
        assert d["recovery"][0] == 1.0, f"{name} must recover at rate 0"
    # screening packs 2 bits/nt vs rotating's ~1.58, at the same ECC
    assert fam["screening+RS"]["bits_per_nt"] > fam["rotating+RS"]["bits_per_nt"]


def test_cost_grid_shape_and_monotonic():
    data = bytes(range(256))
    cov, ov = [1, 4, 16], [0.25, 0.5, 1.0]
    g = benchmark.cost_grid(data, rate=0.01, coverages=cov, overheads=ov, trials=3)
    rec = g["recovery"]
    assert len(rec) == 3 and all(len(r) == 3 for r in rec)
    assert all(0.0 <= v <= 1.0 for r in rec for v in r)
    # more coverage helps (weakly) at the highest redundancy
    assert rec[-1][-1] >= rec[0][-1]


def test_indel_confinement_beats_naive():
    r = benchmark.indel_confinement(bases_len=256, del_rates=[0.0, 0.02],
                                    period=16, trials=4)
    assert r["with_markers"][0] == 1.0 and r["without"][0] == 1.0  # clean at rate 0
    assert r["with_markers"][1] > r["without"][1], "markers must confine deletions"


def test_indel_pipeline_integration_helps():
    r = benchmark.indel_pipeline(del_rates=[0.0, 0.01], coverage=3, trials=4)
    assert r["detect_erase"][0] == 1.0 and r["marker_inner"][0] == 1.0  # clean
    assert r["marker_inner"][1] > r["detect_erase"][1], "integration must help at 1%"


def test_code_rate_keys():
    oligos, meta = ecc.encode(DATA, data_bytes=8, parity_records=4)
    cr = benchmark.code_rate(len(DATA), oligos, meta)
    for k in ("bits_per_nt", "code_rate", "n_oligos", "total_nt"):
        assert k in cr, f"missing code_rate key: {k}"
    assert cr["bits_per_nt"] > 0


def test_constraint_stats_per_oligo():
    oligos, meta = ecc.encode(DATA, data_bytes=8, parity_records=4)
    cs = benchmark.constraint_stats(oligos)
    assert "gc_mean" in cs and "max_homopolymer_run" in cs
    # rotating code guarantees run == 1 per oligo
    assert cs["max_homopolymer_run"] == 1


if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
        print(f"ok  {t.__name__}")
    print(f"\n{len(tests)} checks passed")
