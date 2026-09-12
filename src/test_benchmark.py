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
