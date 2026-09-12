"""Self-checks for the screening codec (M2b). Run: python src/test_screen_codec.py"""
import os
import random
import sys

sys.path.insert(0, os.path.dirname(__file__))
import constraints
import screen_codec

GC = (0.4, 0.6)
MAXRUN = 3


def test_roundtrip_various_lengths():
    rng = random.Random(0)
    for n in (0, 1, 15, 16, 17, 100, 500):
        data = bytes(rng.randrange(256) for _ in range(n))
        oligos, meta = screen_codec.encode(data, gc_range=GC, max_run=MAXRUN)
        assert screen_codec.decode(oligos, meta) == data, f"roundtrip failed at n={n}"


def test_output_satisfies_gc_and_homopolymer():
    rng = random.Random(1)
    data = bytes(rng.randrange(256) for _ in range(500))
    oligos, meta = screen_codec.encode(data, gc_range=GC, max_run=MAXRUN)
    assert oligos
    for o in oligos:
        assert set(o) <= set("ACGT")
        c = constraints.check(o, gc_range=GC, max_run=MAXRUN)
        assert c["gc_ok"], f"GC {c['gc_fraction']:.3f} outside {GC}"
        assert c["homopolymer_ok"], f"run {c['max_homopolymer_run']} exceeds {MAXRUN}"


def test_tight_gc_window_still_roundtrips():
    rng = random.Random(2)
    data = bytes(rng.randrange(256) for _ in range(300))
    tight = (0.45, 0.55)
    oligos, meta = screen_codec.encode(data, gc_range=tight, max_run=MAXRUN)
    for o in oligos:
        assert constraints.check(o, gc_range=tight, max_run=MAXRUN)["gc_ok"]
    assert screen_codec.decode(oligos, meta) == data


def test_density_near_two_bits_per_nt():
    # Large payload per oligo -> per-oligo overhead (seed) amortizes toward 2 bits/nt.
    data = bytes(range(256)) * 4  # 1024 bytes
    oligos, meta = screen_codec.encode(data, payload_bytes=128, gc_range=GC, max_run=MAXRUN)
    nt = sum(len(o) for o in oligos)
    bits_per_nt = len(data) * 8 / nt
    assert 1.7 < bits_per_nt <= 2.0, f"expected ~2 bits/nt, got {bits_per_nt:.3f}"


def test_order_independent_reassembly():
    rng = random.Random(3)
    data = bytes(rng.randrange(256) for _ in range(500))
    oligos, meta = screen_codec.encode(data, gc_range=GC, max_run=MAXRUN)
    shuffled = oligos[:]
    rng.shuffle(shuffled)
    assert screen_codec.decode(shuffled, meta) == data


def test_decode_one_recovers_index():
    data = bytes(range(100))
    oligos, meta = screen_codec.encode(data, payload_bytes=16, gc_range=GC, max_run=MAXRUN)
    idxs = sorted(screen_codec.decode_one(o, meta)[0] for o in oligos)
    assert idxs == list(range(len(oligos)))


if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
        print(f"ok  {t.__name__}")
    print(f"\n{len(tests)} checks passed")
