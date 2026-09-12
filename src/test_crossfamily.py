"""S3 self-checks: every {base codec} x {ECC layer} family works end to end.

Proves the pluggable-base refactor: the RS and fountain ECC layers both accept
either base codec (rotating or screening). Run: python src/test_crossfamily.py
"""
import os
import random
import sys

sys.path.insert(0, os.path.dirname(__file__))
import codec
import ecc
import fountain
import screen_codec

FAMILIES = [
    ("rotating+RS", codec, ecc),
    ("rotating+fountain", codec, fountain),
    ("screening+RS", screen_codec, ecc),
    ("screening+fountain", screen_codec, fountain),
]


def _encode(base, layer, data):
    if layer is ecc:  # ~2x redundancy: parity_records == data records
        return layer.encode(data, data_bytes=16, parity_records=16, block_records=64, base=base)
    return layer.encode(data, data_bytes=16, overhead=2.0, base=base)  # fountain 3x droplets


def test_all_families_roundtrip_clean():
    data = bytes(range(256))
    for name, base, layer in FAMILIES:
        oligos, meta = _encode(base, layer, data)
        assert layer.decode(oligos, meta, base=base) == data, f"{name} clean roundtrip"


def test_all_families_recover_dropout():
    data = bytes(range(256))
    for name, base, layer in FAMILIES:
        rng = random.Random(0)
        oligos, meta = _encode(base, layer, data)
        keep = [o for o in oligos if rng.random() > 0.2]  # lose ~20%
        assert layer.decode(keep, meta, base=base) == data, f"{name} dropout recovery"


if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
        print(f"ok  {t.__name__}")
    print(f"\n{len(tests)} checks passed")
