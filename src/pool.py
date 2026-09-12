"""Random access over a mixed file pool (S5).

Real DNA-storage archives keep many files in one physical pool and retrieve a
single file by amplifying only its oligos with file-specific PCR primers. This
module models that at the codec level: every oligo is prefixed with a short
file-id barcode, so retrieving a file is filtering the pool by barcode and
decoding only that subset, without touching the other files. It reuses the
existing codec/ECC stack (any base codec, RS or fountain), so random access is a
thin addressing layer, not a new codec.

Barcode errors from the channel simply drop an oligo from its file's subset, which
the ECC layer recovers like any other dropout.
"""
import os
import random
import sys

sys.path.insert(0, os.path.dirname(__file__))
import codec
import ecc

BASES = "ACGT"


def _barcode(file_id, barcode_len):
    """A file id as `barcode_len` bases, two bits per base."""
    return "".join(BASES[(file_id >> (2 * (barcode_len - 1 - i))) & 3]
                   for i in range(barcode_len))


def write_pool(files, base=codec, layer=ecc, barcode_len=6, seed=0, **enc_kwargs):
    """Encode a dict of {file_id: bytes} into one shuffled, barcoded oligo pool.

    Returns (pool, pool_meta) where pool_meta carries the barcode width and each
    file's decode meta.
    """
    pool, metas = [], {}
    for file_id, data in files.items():
        oligos, meta = layer.encode(data, base=base, **enc_kwargs)
        metas[file_id] = meta
        bc = _barcode(file_id, barcode_len)
        pool.extend(bc + o for o in oligos)
    random.Random(seed).shuffle(pool)  # mix the files, as a real pool is mixed
    return pool, {"barcode_len": barcode_len, "files": metas}


def select(pool, file_id, pool_meta):
    """The oligos in `pool` that belong to `file_id`, with the barcode stripped."""
    bl = pool_meta["barcode_len"]
    bc = _barcode(file_id, bl)
    return [o[bl:] for o in pool if o[:bl] == bc]


def read_file(pool, file_id, pool_meta, base=codec, layer=ecc):
    """Retrieve one file from the pool: filter by barcode, then decode."""
    return layer.decode(select(pool, file_id, pool_meta),
                        pool_meta["files"][file_id], base=base)


if __name__ == "__main__":
    files = {i: bytes((i * 41 + j) % 256 for j in range(200)) for i in range(1, 6)}
    p, meta = write_pool(files, data_bytes=16, parity_records=8)
    fid = 3
    touched = len(select(p, fid, meta))
    ok = read_file(p, fid, meta) == files[fid]
    print(f"pool     : {len(files)} files, {len(p)} oligos total")
    print(f"read file {fid}: touched {touched} oligos "
          f"({touched / len(p):.0%} of the pool), correct={ok}")
