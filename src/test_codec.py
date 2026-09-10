"""M1 self-checks for the constraint-aware rotating codec.

Assert-based, no framework. Run: python src/test_codec.py
"""
import os
import random
import sys

sys.path.insert(0, os.path.dirname(__file__))

import constraints
import codec


def test_constraints_measure():
    assert constraints.gc_fraction("GCGC") == 1.0
    assert constraints.gc_fraction("ATAT") == 0.0
    assert constraints.gc_fraction("") == 0.0
    assert abs(constraints.gc_fraction("ACGT") - 0.5) < 1e-9
    assert constraints.max_homopolymer_run("AAAC") == 3
    assert constraints.max_homopolymer_run("ACGT") == 1
    assert constraints.max_homopolymer_run("") == 0
    assert constraints.has_forbidden_motif("ACGGATCC", ("GGATCC",)) is True
    assert constraints.has_forbidden_motif("ACGT", ("GGATCC",)) is False


def test_roundtrip_various_lengths():
    rng = random.Random(0)
    for n in (0, 1, 15, 16, 17, 100, 1000):
        data = bytes(rng.randrange(256) for _ in range(n))
        oligos, meta = codec.encode(data)
        assert codec.decode(oligos, meta) == data, f"roundtrip failed at n={n}"


def test_output_is_dna_and_homopolymer_bounded():
    rng = random.Random(1)
    data = bytes(rng.randrange(256) for _ in range(500))
    oligos, meta = codec.encode(data)
    assert oligos, "expected at least one oligo for non-empty data"
    for o in oligos:
        assert set(o) <= set("ACGT"), "oligo has non-ACGT chars"
        # rotating code: consecutive bases always differ -> max run is 1
        assert constraints.max_homopolymer_run(o) == 1


def test_order_independent_reassembly():
    rng = random.Random(2)
    data = bytes(rng.randrange(256) for _ in range(1000))
    oligos, meta = codec.encode(data)
    shuffled = oligos[:]
    rng.shuffle(shuffled)
    assert codec.decode(shuffled, meta) == data


def test_gc_reported_in_range():
    rng = random.Random(3)
    data = bytes(rng.randrange(256) for _ in range(1000))
    oligos, meta = codec.encode(data)
    whole = "".join(oligos)
    gc = constraints.gc_fraction(whole)
    assert 0.0 < gc < 1.0, f"GC out of range: {gc}"


def test_indices_distinct_and_cover_chunks():
    rng = random.Random(4)
    data = bytes(rng.randrange(256) for _ in range(100))
    oligos, meta = codec.encode(data, payload_bytes=16)
    idxs = sorted(codec.oligo_index(o, meta) for o in oligos)
    assert idxs == list(range(len(oligos))), f"indices not 0..k-1: {idxs}"


if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
        print(f"ok  {t.__name__}")
    print(f"\n{len(tests)} checks passed")
