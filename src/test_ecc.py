"""M2 self-checks for the Reed-Solomon error-correction layer.

Assert-based, no framework. Run: python src/test_ecc.py

Corruption here is applied at the DNA level (drop or mutate oligos), which is
what the M3 channel will do properly. M2 uses a detect-and-erase model: an inner
CRC (over index + data) detects any damaged oligo, which the outer striped
Reed-Solomon code then reconstructs. Substitutions and indels both reduce to
"bad oligo -> erasure". M2 validates:
- whole dropped oligos are recovered by the outer code,
- corrupted oligos (payload or index) fail the CRC and are recovered,
- a mix of drops and corruptions up to the parity budget is recovered,
- exceeding the parity budget fails loudly.

In-place correction of small substitutions (saving an erasure slot) needs
error-localized encoding and is a documented later enhancement, not M2.
"""
import os
import random
import sys

sys.path.insert(0, os.path.dirname(__file__))

import codec
import ecc

CFG = dict(data_bytes=8, parity_records=4, block_records=64)


def _garble(oligo):
    """Replace an oligo with an unrelated valid rotating string of equal length,
    enough byte damage that the inner code flags it as an erasure."""
    rng = random.Random(hash(oligo) & 0xFFFF)
    out, prev = [], "A"
    for _ in range(len(oligo)):
        b = rng.choice([x for x in "ACGT" if x != prev])
        out.append(b)
        prev = b
    return "".join(out)


def test_roundtrip_clean():
    rng = random.Random(0)
    for n in (0, 1, 8, 40, 200, 1000):
        data = bytes(rng.randrange(256) for _ in range(n))
        oligos, meta = ecc.encode(data, **CFG)
        assert ecc.decode(oligos, meta) == data, f"clean roundtrip failed at n={n}"


def test_decode_records_reports_gaps():
    data = bytes(range(40))
    oligos, meta = ecc.encode(data, **CFG)
    survivors = codec.decode_records(oligos[:-2], meta)  # drop last two
    assert isinstance(survivors, dict)
    dropped = {len(oligos) - 2, len(oligos) - 1}
    assert dropped.isdisjoint(survivors.keys()), "dropped indices must be absent"
    assert set(survivors) == set(range(len(oligos))) - dropped


def test_recovers_dropped_oligos():
    data = bytes(range(40))  # 5 data records, +4 parity = 9 oligos, one block
    oligos, meta = ecc.encode(data, **CFG)
    rng = random.Random(1)
    keep = oligos[:]
    for i in rng.sample(range(len(oligos)), CFG["parity_records"]):  # drop exactly P
        keep[i] = None
    keep = [o for o in keep if o is not None]
    assert ecc.decode(keep, meta) == data


def test_recovers_corrupted_oligos():
    data = bytes(range(40))
    oligos, meta = ecc.encode(data, **CFG)
    rng = random.Random(2)
    idx = rng.sample(range(len(oligos)), CFG["parity_records"])  # garble P of them
    corrupted = [_garble(o) if i in idx else o for i, o in enumerate(oligos)]
    assert ecc.decode(corrupted, meta) == data


def test_recovers_mixed_drops_and_corruptions():
    data = bytes(range(40))  # 5 data + 4 parity = 9 oligos, one block
    oligos, meta = ecc.encode(data, **CFG)
    rng = random.Random(5)
    victims = rng.sample(range(len(oligos)), CFG["parity_records"])  # total damage == P
    drop, garble = set(victims[:2]), set(victims[2:])
    mixed = [
        _garble(o) if i in garble else o
        for i, o in enumerate(oligos)
        if i not in drop
    ]
    assert ecc.decode(mixed, meta) == data


def test_exceeds_capacity_raises():
    data = bytes(range(40))
    oligos, meta = ecc.encode(data, **CFG)
    keep = [o for i, o in enumerate(oligos) if i >= CFG["parity_records"] + 1]  # drop P+1
    try:
        ecc.decode(keep, meta)
        assert False, "expected failure when erasures exceed parity budget"
    except ecc.RecoveryError:
        pass


if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
        print(f"ok  {t.__name__}")
    print(f"\n{len(tests)} checks passed")
