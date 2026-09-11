"""M3 self-checks for the error channel. Run: python src/test_channel.py"""
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
import channel

OLIGO = "ACGTACGTAC"
POOL = ["ACGTAC", "GTACGT", "TACGTA", "CGTACG"]


def test_zero_rates_identity():
    out = channel.corrupt(POOL, seed=0)
    assert sorted(out) == sorted(POOL), "zero-rate channel must return oligos unchanged"


def test_substitution_preserves_length_changes_bases():
    (read,) = channel.corrupt([OLIGO], p_sub=1.0, seed=0)
    assert len(read) == len(OLIGO), "substitution must not change length"
    assert all(r != o for r, o in zip(read, OLIGO)), "p_sub=1 must change every base"
    assert set(read) <= set("ACGT")


def test_deletion_empties():
    (read,) = channel.corrupt([OLIGO], p_del=1.0, seed=0)
    assert read == "", "p_del=1 must delete every base"


def test_insertion_lengthens():
    (read,) = channel.corrupt([OLIGO], p_ins=1.0, seed=0)
    assert len(read) == 2 * len(OLIGO), "p_ins=1 inserts one base per position"
    assert set(read) <= set("ACGT")


def test_dropout_drops_all():
    assert channel.corrupt(POOL, p_drop=1.0, seed=0) == []


def test_coverage_makes_multiple_reads():
    out = channel.corrupt([OLIGO], coverage=3, seed=0)
    assert out == [OLIGO, OLIGO, OLIGO], "coverage=3 at zero rate -> 3 clean copies"


def test_determinism():
    a = channel.corrupt(POOL, p_sub=0.1, p_ins=0.05, p_del=0.05, seed=7)
    b = channel.corrupt(POOL, p_sub=0.1, p_ins=0.05, p_del=0.05, seed=7)
    assert a == b, "same seed must give identical output"


if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
        print(f"ok  {t.__name__}")
    print(f"\n{len(tests)} checks passed")
