"""M2 error-correction layer for the DNA-storage codec.

Detect-and-erase model. Every record carries an inner CRC over (index + data);
any damaged oligo (substitution, garble, corrupted index) fails the CRC and is
treated as an erasure. An outer striped Reed-Solomon code (systematic, GF(256))
adds parity records so up to `parity_records` erased or damaged oligos per block
are reconstructed. A dropped oligo is a missing index, also an erasure. This is
how DNA-storage codecs survive substitutions and, at M3, insertions/deletions
(a de-synchronized read fails its CRC and is recovered as an erasure).

Wraps the M1 codec's `bytes <-> oligos` interface without changing its encoding.

ponytail: no in-place inner correction of small substitutions; every damaged
oligo spends one erasure slot. Add error-localized encoding + inner RS correction
if the recovery-vs-redundancy curve needs it.
"""
import os
import sys
import zlib

sys.path.insert(0, os.path.dirname(__file__))
import codec
import markercode

from reedsolo import RSCodec, ReedSolomonError

CRC_BYTES = 4
RS_MAX = 255  # GF(256): data records + parity records per block must fit


class RecoveryError(Exception):
    """Raised when erasures in a block exceed the parity budget."""


def _crc(index, data):
    return zlib.crc32(index.to_bytes(4, "big") + data) & 0xFFFFFFFF


def encode(data, data_bytes=8, parity_records=4, block_records=64, base=codec,
           inner_nsym=0, marker_period=0):
    """Encode bytes into error-corrected DNA oligos plus out-of-band meta.

    `base` is the base codec module that maps records to oligos (rotating `codec`
    or `screen_codec`); it must expose `encode(data, payload_bytes)` and
    `decode_one(oligo, meta)`. Detection is a fixed 4-byte CRC per record.

    Two optional inner layers make the codec survive insertions/deletions instead
    of erasing the whole oligo (both default off, keeping the plain detect-and-erase
    behavior):
      - `inner_nsym` > 0 adds that many Reed-Solomon parity bytes per oligo, which
        correct the residual substitution-like errors a resynced indel leaves behind;
      - `marker_period` > 0 wraps each oligo with resync markers (see `markercode`)
        so an indel's damage is confined to one run rather than desynchronizing the
        whole oligo.
    Used together they form a concatenated code: markers turn an indel into a few
    local byte errors, and the inner RS corrects them, so the oligo is recovered
    rather than lost.

    This requires an error-LOCALIZED base codec, i.e. `base=screen_codec` (2 bits/nt,
    per-byte packing). The rotating `codec` encodes each record as one base-3 big
    integer, so a single base error propagates through base-256 carries into ~half the
    record's bytes, which the inner RS cannot repair; with the rotating codec these
    layers do not help. This mirrors the M2 note on why in-place correction needs
    localized coding.
    """
    D, P = data_bytes, parity_records
    if block_records + P > RS_MAX:
        raise ValueError(f"block_records + parity_records must be <= {RS_MAX}")

    total_len = len(data)
    padded = data + b"\x00" * ((-len(data)) % D) if data else b""
    data_recs = [padded[i:i + D] for i in range(0, len(padded), D)]
    blocks = [data_recs[i:i + block_records] for i in range(0, len(data_recs), block_records)]

    outer = RSCodec(P)
    records = []          # ordered: block0 data, block0 parity, block1 data, ...
    block_sizes = []
    for blk in blocks:
        K = len(blk)
        block_sizes.append(K)
        parity = [bytearray(D) for _ in range(P)]
        for j in range(D):                       # striped RS, one codeword per column
            enc = outer.encode(bytes(rec[j] for rec in blk))  # length K + P, systematic
            for p in range(P):
                parity[p][j] = enc[K + p]
        records.extend(blk)
        records.extend(bytes(pr) for pr in parity)

    # inner CRC over (global index + record), then optional inner RS parity
    inner = RSCodec(inner_nsym) if inner_nsym else None
    framed = []
    for gi, rec in enumerate(records):
        f = bytes(rec) + _crc(gi, bytes(rec)).to_bytes(CRC_BYTES, "big")
        framed.append(bytes(inner.encode(f)) if inner else f)

    record_len = D + CRC_BYTES + inner_nsym
    blob = b"".join(framed)
    oligos, cmeta = base.encode(blob, payload_bytes=record_len)
    oligo_len = len(oligos[0]) if oligos else 0
    if marker_period:
        oligos = [markercode.encode(o, marker_period) for o in oligos]
    meta = {
        **cmeta,
        "ecc": {
            "data_bytes": D,
            "parity_records": P,
            "block_sizes": block_sizes,
            "total_len": total_len,
            "record_len": record_len,
            "inner_nsym": inner_nsym,
            "marker_period": marker_period,
            "oligo_len": oligo_len,
        },
    }
    return oligos, meta


def decode(oligos, meta, base=codec):
    """Recover the original bytes from oligos (any order, with erasures).

    `base` must be the same base codec used to encode.
    """
    e = meta["ecc"]
    D, P = e["data_bytes"], e["parity_records"]
    inner_nsym, marker_period = e.get("inner_nsym", 0), e.get("marker_period", 0)
    oligo_len = e.get("oligo_len", 0)
    inner = RSCodec(inner_nsym) if inner_nsym else None
    outer = RSCodec(P)

    # Decode reads individually and keep the first CRC-valid read per index.
    # With sequencing coverage (multiple reads per oligo), one clean read is
    # enough to fill a slot, so a damaged read of an otherwise-covered oligo
    # does not cost an erasure. When enabled, markers resync each read (confining
    # indels) and the inner RS corrects the residual byte errors before the CRC.
    good = {}
    for oligo in oligos:
        if marker_period:
            oligo = markercode.decode(oligo, oligo_len, marker_period)
        idx, framed = base.decode_one(oligo, meta)
        if idx in good:
            continue
        if inner:
            try:
                framed = bytes(inner.decode(framed)[0])
            except ReedSolomonError:
                continue  # inner code could not repair this read
        rec, crc = framed[:D], framed[D:D + CRC_BYTES]
        if _crc(idx, rec) == int.from_bytes(crc, "big"):
            good[idx] = rec

    out = bytearray()
    gi = 0
    for K in e["block_sizes"]:
        N = K + P
        present, erase = {}, []
        for slot in range(N):
            g = gi + slot
            if g in good:
                present[slot] = good[g]
            else:
                erase.append(slot)
        if len(erase) > P:
            raise RecoveryError(
                f"block at index {gi}: {len(erase)} erasures exceed parity budget {P}"
            )
        data_recs = [bytearray(D) for _ in range(K)]
        for j in range(D):
            word = bytearray(N)
            for slot, rec in present.items():
                word[slot] = rec[j]
            msg = outer.decode(bytes(word), erase_pos=erase)[0] if erase \
                else outer.decode(bytes(word))[0]
            for i in range(K):
                data_recs[i][j] = msg[i]
        for rec in data_recs:
            out += rec
        gi += N
    return bytes(out[: e["total_len"]])


if __name__ == "__main__":
    import constraints

    msg = b"DNA storage with Reed-Solomon error correction!"
    oligos, meta = encode(msg, data_bytes=8, parity_records=4)
    # simulate a dropout and a garble within budget
    survivors = oligos[:-1]                      # drop one
    survivors[0] = "".join(reversed(survivors[0]))  # corrupt one
    back = decode(survivors, meta)
    whole = "".join(oligos)
    print(f"message : {msg!r}")
    print(f"oligos  : {len(oligos)} (payload+parity), {len(whole)} nt")
    print(f"decoded : {back!r}  (match={back == msg})")
    print(f"overhead: {len(oligos)}/{len(oligos) - meta['ecc']['parity_records']} records "
          f"(+{meta['ecc']['parity_records']} parity per block)")
    # homopolymer runs are a per-oligo property (each oligo is one molecule);
    # measuring over the joined pool would invent runs at oligo boundaries.
    per_oligo_run = max(constraints.max_homopolymer_run(o) for o in oligos)
    print(f"GC      : {constraints.gc_fraction(whole):.3f}  max run (per oligo): {per_oligo_run}")
