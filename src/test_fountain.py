"""Self-checks for the fountain (LT) codec (S2). Run: python src/test_fountain.py"""
import os
import random
import sys

sys.path.insert(0, os.path.dirname(__file__))
import fountain

CFG = dict(data_bytes=16, overhead=3.0)  # generous overhead -> deterministic peeling


def test_roundtrip_clean():
    rng = random.Random(0)
    for n in (0, 1, 16, 100, 500):
        data = bytes(rng.randrange(256) for _ in range(n))
        oligos, meta = fountain.encode(data, **CFG)
        assert fountain.decode(oligos, meta) == data, f"clean roundtrip failed at n={n}"


def test_more_droplets_than_segments():
    data = bytes(range(256))  # 16 segments at data_bytes=16
    oligos, meta = fountain.encode(data, **CFG)
    assert len(oligos) > meta["k"], "fountain must emit more droplets than segments"


def test_recovers_from_dropout():
    data = bytes(range(256))
    oligos, meta = fountain.encode(data, **CFG)
    rng = random.Random(1)
    # drop a third of the oligos; enough droplets remain to peel
    keep = [o for o in oligos if rng.random() > 0.33]
    assert fountain.decode(keep, meta) == data


def test_recovers_from_corruption():
    data = bytes(range(256))
    oligos, meta = fountain.encode(data, **CFG)
    rng = random.Random(2)

    def garble(o):
        out, prev = [], "A"
        for _ in range(len(o)):
            b = rng.choice([x for x in "ACGT" if x != prev]); out.append(b); prev = b
        return "".join(out)

    corrupted = [garble(o) if rng.random() < 0.2 else o for o in oligos]
    assert fountain.decode(corrupted, meta) == data


def test_too_few_droplets_raises():
    data = bytes(range(256))
    oligos, meta = fountain.encode(data, **CFG)
    survivors = oligos[: meta["k"] - 2]  # fewer than K -> unrecoverable
    try:
        fountain.decode(survivors, meta)
        assert False, "expected RecoveryError with fewer than K droplets"
    except fountain.RecoveryError:
        pass


if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
        print(f"ok  {t.__name__}")
    print(f"\n{len(tests)} checks passed")
