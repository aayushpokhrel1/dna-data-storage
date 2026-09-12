"""Self-checks for random access over a mixed file pool (S5).
Run: python src/test_pool.py
"""
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
import channel
import pool


def _files():
    return {
        1: b"the quick brown fox " * 12,
        2: bytes(range(200)),
        3: b"third file, different content entirely " * 6,
        4: bytes((i * 37 + 5) % 256 for i in range(150)),
    }


def test_write_read_roundtrip_clean():
    files = _files()
    p, meta = pool.write_pool(files, data_bytes=16, parity_records=8)
    for fid, data in files.items():
        assert pool.read_file(p, fid, meta) == data, f"file {fid} roundtrip"


def test_random_access_touches_only_target():
    files = _files()
    p, meta = pool.write_pool(files, data_bytes=16, parity_records=8)
    sel = pool.select(p, 2, meta)
    assert 0 < len(sel) < len(p), "random access must read a strict subset of the pool"
    # the touched share is roughly one file out of four
    assert len(sel) / len(p) < 0.5


def test_cross_file_isolation():
    files = _files()
    p, meta = pool.write_pool(files, data_bytes=16, parity_records=8)
    assert pool.read_file(p, 1, meta) == files[1]
    assert pool.read_file(p, 1, meta) != files[2]


def test_recovers_under_channel():
    files = _files()
    p, meta = pool.write_pool(files, data_bytes=16, parity_records=12)
    reads = channel.corrupt(p, p_sub=0.005, p_drop=0.1, coverage=6, seed=0)
    for fid, data in files.items():
        assert pool.read_file(reads, fid, meta) == data, f"file {fid} under channel"


if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
        print(f"ok  {t.__name__}")
    print(f"\n{len(tests)} checks passed")
