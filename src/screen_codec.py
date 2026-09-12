"""Screening codec (M2b): higher-density constraint satisfaction by rejection.

DNA Fountain-style (Erlich & Zielinski 2017). Bits pack two per base
(A=00, C=01, G=10, T=11), so the payload region is ~2 bits/nt, above the rotating
codec's ~1.58. A per-chunk seed drives a hash keystream that scrambles the record
`[index + payload]`; the encoder tries seeds until the packed oligo satisfies BOTH
a GC window and the homopolymer bound, then stores the seed in the clear. So this
codec guarantees GC in `gc_range` AND max run <= `max_run` BY CONSTRUCTION, the GC
guarantee the rotating codec lacks. The cost is compute (rejections), not
nucleotides.

Same interface as the rotating codec (`encode`/`decode`/`decode_one`, indexed
oligos, ECC-agnostic), so it drops into the same channel and benchmark.

meta is out-of-band, as in the rotating codec.
"""
import hashlib

BASES = "ACGT"
_VAL = {b: i for i, b in enumerate(BASES)}


def _pack(data):
    """Bytes -> bases, two bits per base, most-significant pair first."""
    out = []
    for byte in data:
        out.append(BASES[(byte >> 6) & 3])
        out.append(BASES[(byte >> 4) & 3])
        out.append(BASES[(byte >> 2) & 3])
        out.append(BASES[byte & 3])
    return "".join(out)


def _unpack(seq):
    """Bases -> bytes (inverse of _pack; len(seq) must be a multiple of 4)."""
    out = bytearray()
    for i in range(0, len(seq), 4):
        b = (_VAL[seq[i]] << 6) | (_VAL[seq[i + 1]] << 4) \
            | (_VAL[seq[i + 2]] << 2) | _VAL[seq[i + 3]]
        out.append(b)
    return bytes(out)


def _keystream(seed, n):
    out = bytearray()
    ctr = 0
    while len(out) < n:
        out += hashlib.sha256(seed.to_bytes(8, "big") + ctr.to_bytes(8, "big")).digest()
        ctr += 1
    return bytes(out[:n])


def _scramble(data, seed):
    """XOR data with a seed-derived keystream. Its own inverse."""
    return bytes(b ^ k for b, k in zip(data, _keystream(seed, len(data))))


def _run_ok(seq, max_run):
    run = 1
    for i in range(1, len(seq)):
        run = run + 1 if seq[i] == seq[i - 1] else 1
        if run > max_run:
            return False
    return True


def _gc_ok(seq, gc_range):
    gc = sum(b in "GC" for b in seq) / len(seq) if seq else 0.0
    return gc_range[0] <= gc <= gc_range[1]


class ScreenError(Exception):
    """No seed within max_seed produced a constraint-satisfying oligo."""


def encode(data, payload_bytes=16, index_bytes=3, seed_bases=8,
           gc_range=(0.4, 0.6), max_run=3, max_seed=1 << 16):
    """Encode bytes into constraint-satisfying DNA oligos plus out-of-band meta."""
    meta = {
        "payload_bytes": payload_bytes,
        "index_bytes": index_bytes,
        "seed_bases": seed_bases,
        "total_len": len(data),
    }
    oligos = []
    for idx, start in enumerate(range(0, len(data), payload_bytes)):
        chunk = data[start:start + payload_bytes]
        chunk = chunk + b"\x00" * (payload_bytes - len(chunk))  # zero-pad last
        record = idx.to_bytes(index_bytes, "big") + chunk
        for seed in range(max_seed):
            oligo = _pack(_seed_bytes(seed, seed_bases)) + _pack(_scramble(record, seed))
            if _gc_ok(oligo, gc_range) and _run_ok(oligo, max_run):
                oligos.append(oligo)
                break
        else:
            raise ScreenError(f"no valid seed < {max_seed} for chunk {idx}")
    return oligos, meta


def _seed_bytes(seed, seed_bases):
    """Seed value as the byte string that _pack turns into exactly seed_bases bases."""
    return seed.to_bytes(seed_bases // 4, "big")


def decode_one(oligo, meta):
    """Decode one oligo to (index, payload_bytes).

    Length-tolerant: a channel read may be shortened or lengthened by indels, so
    the seed region is padded to width and the body truncated to a whole number of
    bases. A malformed read still yields a record, which the caller's CRC rejects.
    """
    sb = meta["seed_bases"]
    seed_region = (oligo[:sb] + "A" * sb)[:sb]
    body = oligo[sb:]
    body = body[: (len(body) // 4) * 4]
    seed = int.from_bytes(_unpack(seed_region), "big")
    record = _scramble(_unpack(body), seed)
    ib = meta["index_bytes"]
    idx = int.from_bytes(record[:ib], "big") if len(record) >= ib else 0
    return idx, record[ib:ib + meta["payload_bytes"]]


def decode(oligos, meta):
    """Reassemble the original bytes from oligos (any order)."""
    records = dict(decode_one(o, meta) for o in oligos)
    if not records:
        return b""
    data = b"".join(records[i] for i in range(len(records)))
    return data[: meta["total_len"]]


if __name__ == "__main__":
    import constraints

    msg = b"Screening codec: 2 bits/nt with a guaranteed GC window."
    oligos, meta = encode(msg, payload_bytes=16)
    back = decode(oligos, meta)
    whole = "".join(oligos)
    print(f"message : {msg!r}")
    print(f"decoded : {back!r}  (match={back == msg})")
    print(f"density : {len(msg) * 8 / len(whole):.3f} bits/nt")
    gcs = [constraints.gc_fraction(o) for o in oligos]
    runs = max(constraints.max_homopolymer_run(o) for o in oligos)
    print(f"GC      : {min(gcs):.3f}-{max(gcs):.3f} (per oligo)   max run: {runs}")
