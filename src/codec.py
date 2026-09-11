"""Constraint-aware rotating-code codec for DNA data storage (M1).

Maps bytes to DNA oligos and back. The rotating rule picks each base from the
three that differ from the previous base, so every adjacent pair differs and the
maximum homopolymer run is 1 by construction (~1.58 bits/nt). Payload is chunked
into indexed oligos for an addressable pool; decode reassembles by index, so
oligos may arrive shuffled or partial.

ECC-agnostic: the interface is bytes <-> oligos. The error-correction layer (M2)
wraps this without changing it.

meta is passed out-of-band for M1.
ponytail: no self-describing in-DNA header yet; add one when the pool must be
parameter-free (random access from a bare oligo set).
"""
BASES = "ACGT"
SEED_BASE = "A"  # virtual "previous" base defining the first real base


def _alts(prev):
    """The three bases that differ from prev, in fixed BASES order."""
    return [b for b in BASES if b != prev]


def _trits_to_dna(trits, prev=SEED_BASE):
    out = []
    for t in trits:
        b = _alts(prev)[t]
        out.append(b)
        prev = b
    return "".join(out)


def _dna_to_trits(seq, prev=SEED_BASE):
    trits = []
    for b in seq:
        alts = _alts(prev)
        # Corruption-tolerant: a substitution can make a base equal its
        # predecessor (violating the rotating invariant) or otherwise land
        # outside the alternatives. Emit a placeholder trit and keep going;
        # the resulting byte error is caught by the M2 inner code.
        trits.append(alts.index(b) if b in alts else 0)
        prev = b
    return trits


def _trit_width(n_bytes):
    """Smallest T with 3**T >= 256**n_bytes (exact integer check)."""
    target = 256 ** n_bytes
    t, val = 0, 1
    while val < target:
        val *= 3
        t += 1
    return t


def _int_to_trits(n, width):
    trits = []
    for _ in range(width):
        n, r = divmod(n, 3)
        trits.append(r)
    assert n == 0, "value does not fit in the given trit width"
    return trits[::-1]


def _trits_to_int(trits):
    n = 0
    for t in trits:
        n = n * 3 + t
    return n


def encode(data, payload_bytes=16, index_trits=20, seed=SEED_BASE):
    """Encode bytes into a list of DNA oligos plus out-of-band meta."""
    pt_width = _trit_width(payload_bytes)
    meta = {
        "total_len": len(data),
        "payload_bytes": payload_bytes,
        "index_trits": index_trits,
        "payload_trits": pt_width,
        "seed": seed,
    }
    oligos = []
    for idx, start in enumerate(range(0, len(data), payload_bytes)):
        chunk = data[start:start + payload_bytes]
        chunk = chunk + b"\x00" * (payload_bytes - len(chunk))  # zero-pad last
        trits = (
            _int_to_trits(idx, index_trits)
            + _int_to_trits(int.from_bytes(chunk, "big"), pt_width)
        )
        oligos.append(_trits_to_dna(trits, seed))
    return oligos, meta


def oligo_index(oligo, meta):
    """Recover just the chunk index from one oligo."""
    trits = _dna_to_trits(oligo, meta["seed"])
    return _trits_to_int(trits[: meta["index_trits"]])


def decode_one(oligo, meta):
    """Decode one oligo (one read) to (index, payload_bytes).

    A payload integer wider than payload_bytes (corruption) is masked to fit
    rather than raising, so a bad read still yields a record the caller can
    reject via its own integrity check.
    """
    span = 256 ** meta["payload_bytes"]
    trits = _dna_to_trits(oligo, meta["seed"])
    idx = _trits_to_int(trits[: meta["index_trits"]])
    val = _trits_to_int(trits[meta["index_trits"]:]) % span
    return idx, val.to_bytes(meta["payload_bytes"], "big")


def decode_records(oligos, meta):
    """Recover records keyed by index: {index: payload_bytes}.

    Only present indices appear, so callers (the M2 ECC layer) can see which
    are missing and treat them as erasures. Duplicate reads of one index (from
    sequencing coverage) collapse here (last wins); the ECC layer decodes reads
    individually so it can keep an integrity-valid read instead.
    """
    return dict(decode_one(o, meta) for o in oligos)


def decode(oligos, meta):
    """Reassemble the original bytes from oligos (any order). No error
    correction: assumes every oligo is present and clean."""
    records = decode_records(oligos, meta)
    if not records:
        return b""
    data = b"".join(records[i] for i in range(len(records)))
    return data[: meta["total_len"]]


if __name__ == "__main__":
    import constraints

    msg = b"Hello, DNA data storage!"
    oligos, meta = encode(msg, payload_bytes=8)
    back = decode(oligos, meta)
    whole = "".join(oligos)
    print(f"message  : {msg!r}")
    print(f"oligos   : {oligos}")
    print(f"decoded  : {back!r}  (match={back == msg})")
    print(f"nt total : {len(whole)}  bits/nt={len(msg) * 8 / len(whole):.3f}")
    print(f"GC       : {constraints.gc_fraction(whole):.3f}")
    print(f"max run  : {constraints.max_homopolymer_run(whole)}")
