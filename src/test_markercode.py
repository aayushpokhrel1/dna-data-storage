"""Self-checks for the marker resync inner code (S6).
Run: python src/test_markercode.py
"""
import os
import random
import sys

sys.path.insert(0, os.path.dirname(__file__))
import markercode


def _rotating_bases(n, seed=0):
    """A homopolymer-free base string (like a rotating-codec oligo)."""
    rng = random.Random(seed)
    out, prev = [], "A"
    for _ in range(n):
        b = rng.choice([x for x in "ACGT" if x != prev])
        out.append(b)
        prev = b
    return "".join(out)


def _correct_fraction(a, b):
    """Fraction of positions that match, over the shorter length."""
    n = min(len(a), len(b))
    return sum(x == y for x, y in zip(a[:n], b[:n])) / n if n else 0.0


def test_roundtrip_no_channel():
    for n in (1, 16, 17, 100, 257):
        bases = _rotating_bases(n, seed=n)
        enc = markercode.encode(bases, period=16)
        assert markercode.decode(enc, n, period=16) == bases, f"clean roundtrip n={n}"


def test_single_deletion_is_confined():
    bases = _rotating_bases(256, seed=1)
    enc = markercode.encode(bases, period=16)
    # delete one base in the middle of the encoded stream
    p = len(enc) // 2
    corrupted = enc[:p] + enc[p + 1:]
    resynced = markercode.decode(corrupted, len(bases), period=16)
    # marker resync keeps almost all bases; damage is confined to one run
    assert _correct_fraction(resynced, bases) > 0.9, "resync must confine the deletion"


def test_beats_naive_baseline():
    bases = _rotating_bases(256, seed=2)
    enc = markercode.encode(bases, period=16)
    p = len(enc) // 3
    corrupted = enc[:p] + enc[p + 1:]  # one deletion
    resynced = markercode.decode(corrupted, len(bases), period=16)
    # naive: no markers, a deletion shifts every downstream base
    naive = (bases[:p] + bases[p + 1:])  # same deletion applied to raw data
    assert _correct_fraction(resynced, bases) > _correct_fraction(naive, bases) + 0.2, \
        "marker resync must beat the naive (unmarked) baseline"


if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
        print(f"ok  {t.__name__}")
    print(f"\n{len(tests)} checks passed")
