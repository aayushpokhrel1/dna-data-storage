"""Fountain (LT) codec (S2): a rateless second ECC family for the comparison.

DNA Fountain-style (Erlich & Zielinski 2017). Data is split into K segments. Each
droplet XORs a random subset of segments, with the subset drawn from a Robust
Soliton degree distribution seeded by the droplet number. The droplet number is the
codec index, so no separate seed field is stored. N > K droplets are emitted; a
per-droplet CRC (over index + payload) discards corrupted droplets, so dropouts and
corruptions are both just missing droplets, which is exactly what a fountain code
absorbs. The decoder recomputes each surviving droplet's segment set from its seed
and peels: a degree-1 droplet reveals a segment, which is XORed out of every droplet
that references it, until all K segments resolve.

Same interface as `ecc.py` (`encode(data, ...) -> (oligos, meta)`,
`decode(oligos, meta) -> bytes`, raises `RecoveryError`), so the cross-family
benchmark (S3) can swap the two ECC layers. Wraps the rotating codec for base
mapping, mirroring `ecc.py`.

ponytail: peeling is the naive O(N * K^2) scan; swap in a ripple queue if the
benchmark throughput needs it.
"""
import math
import os
import random
import sys
import zlib

sys.path.insert(0, os.path.dirname(__file__))
import codec

CRC_BYTES = 4


class RecoveryError(Exception):
    """Peeling stalled before all K segments were recovered."""


def _crc(index, data):
    return zlib.crc32(index.to_bytes(4, "big") + data) & 0xFFFFFFFF


def _rsd_cumulative(K, c, delta):
    """Cumulative Robust Soliton distribution over degrees 1..K."""
    rho = [0.0] * (K + 1)
    rho[1] = 1.0 / K
    for d in range(2, K + 1):
        rho[d] = 1.0 / (d * (d - 1))
    tau = [0.0] * (K + 1)
    R = c * math.log(K / delta) * math.sqrt(K) if K > 1 else 0.0
    kr = min(int(K / R), K) if R > 0 else 0  # clamp: K/R can exceed K when R < 1
    for d in range(1, kr):
        tau[d] = R / (d * K)
    if 1 <= kr <= K:
        tau[kr] = R * math.log(R / delta) / K
    w = [rho[d] + tau[d] for d in range(K + 1)]
    Z = sum(w)
    cum, run = [], 0.0
    for d in range(K + 1):
        run += w[d] / Z
        cum.append(run)
    return cum  # cum[d] = P(degree <= d)


def _segments_for(seed, K, cum):
    """Reproduce a droplet's segment set from its seed (encoder and decoder agree)."""
    rng = random.Random(seed)
    r = rng.random()
    degree = next(d for d in range(1, K + 1) if cum[d] >= r)
    return rng.sample(range(K), degree)


def encode(data, data_bytes=16, overhead=0.5, c=0.1, delta=0.5,
           index_trits=20, seed=codec.SEED_BASE):
    """Encode bytes into fountain droplets (DNA oligos) plus out-of-band meta."""
    D = data_bytes
    total_len = len(data)
    padded = data + b"\x00" * ((-len(data)) % D) if data else b""
    segments = [padded[i:i + D] for i in range(0, len(padded), D)]
    K = len(segments)
    meta = {
        "k": K, "data_bytes": D, "c": c, "delta": delta,
        "total_len": total_len, "record_len": D + CRC_BYTES,
    }
    if K == 0:
        oligos, cmeta = codec.encode(b"", payload_bytes=D + CRC_BYTES,
                                     index_trits=index_trits, seed=seed)
        return oligos, {**cmeta, **meta}

    n_droplets = K + max(1, math.ceil(K * overhead))
    cum = _rsd_cumulative(K, c, delta)
    records = []
    for i in range(n_droplets):
        payload = bytearray(D)
        for s in _segments_for(i, K, cum):
            seg = segments[s]
            for j in range(D):
                payload[j] ^= seg[j]
        records.append(bytes(payload) + _crc(i, bytes(payload)).to_bytes(CRC_BYTES, "big"))

    blob = b"".join(records)
    oligos, cmeta = codec.encode(blob, payload_bytes=D + CRC_BYTES,
                                 index_trits=index_trits, seed=seed)
    return oligos, {**cmeta, **meta, "n_droplets": n_droplets}


def decode(oligos, meta):
    """Recover the original bytes by peeling the surviving droplets."""
    D, K = meta["data_bytes"], meta["k"]
    if K == 0:
        return b""
    cum = _rsd_cumulative(K, meta["c"], meta["delta"])

    # Collect CRC-valid droplets as [segment set, payload]; dedupe by seed.
    droplets, seen = [], set()
    for oligo in oligos:
        i, framed = codec.decode_one(oligo, meta)
        if i in seen or i >= meta.get("n_droplets", i + 1):
            continue
        payload, crc = framed[:D], framed[D:D + CRC_BYTES]
        if _crc(i, payload) != int.from_bytes(crc, "big"):
            continue
        seen.add(i)
        droplets.append([set(_segments_for(i, K, cum)), bytearray(payload)])

    solved = [None] * K
    progress = True
    while progress and any(x is None for x in solved):
        progress = False
        for segs, payload in droplets:
            for s in [s for s in segs if solved[s] is not None]:  # peel known segments
                known = solved[s]
                for j in range(D):
                    payload[j] ^= known[j]
                segs.discard(s)
            if len(segs) == 1:
                s = next(iter(segs))
                if solved[s] is None:
                    solved[s] = bytes(payload)
                    segs.discard(s)
                    progress = True

    if any(x is None for x in solved):
        raise RecoveryError(f"peeling stalled: {solved.count(None)}/{K} segments unresolved")
    return b"".join(solved)[: meta["total_len"]]


if __name__ == "__main__":
    # LT peeling is asymptotic; small K needs proportionally more overhead. Use a
    # realistic payload so a modest overhead reliably recovers after dropout.
    msg = (b"Fountain codes give rateless erasure coding: emit as many droplets "
           b"as you like, and any large-enough subset reconstructs the data. ") * 8
    oligos, meta = encode(msg, data_bytes=16, overhead=1.0)
    kept = oligos[: int(len(oligos) * 0.8)]  # lose 20% of droplets
    back = decode(kept, meta)
    print(f"message  : {len(msg)} bytes, {meta['k']} segments")
    print(f"droplets : {meta['n_droplets']} ({meta['n_droplets'] / meta['k']:.2f}x), "
          f"kept {len(kept)} after 20% loss")
    print(f"decoded  : match={back == msg}")
