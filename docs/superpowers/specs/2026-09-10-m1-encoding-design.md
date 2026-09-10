# M1 - Constraint-aware codec (design)

Date: 2026-09-10. Status: approved, implementing.

## Goal

A runnable, tested bits-to-DNA codec that round-trips exactly with no channel and
whose output meets the homopolymer constraint by construction. This is the spine
the M2 error-correction layer wraps and the M3 channel corrupts. Kept ECC-agnostic:
the interface is `bytes <-> DNA oligos`, nothing about coding lives here.

## Decisions (from the user, 2026-09-10)

- **Encoding scheme:** rotating code as the M1 spine (Goldman-style); a higher-density
  screening codec is added later as a second codec and the benchmark compares them.
- **ECC + indel strategy:** deferred to M2, decided with the user. M1 stays ECC-agnostic.
- **Random access:** addressable pool. Payload is chunked into indexed oligos;
  `decode` reassembles by index, so oligos may arrive shuffled or partial.
- **Error rates:** cited from DNA-storage literature (M3), not needed for M1.

## Modules

- `src/constraints.py` - pure predicates and measurement, shared by both codecs and
  the benchmark. No encoding logic.
  - `gc_fraction(seq) -> float`
  - `max_homopolymer_run(seq) -> int`
  - `has_forbidden_motif(seq, motifs=()) -> bool` (stub; wired in when the error model
    names motifs)
  - `check(seq, gc_range, max_run, motifs) -> dict` (measured values + pass flags)

- `src/codec.py` - the rotating-code encoder/decoder.
  - Rotating rule: the next base is chosen from the 3 that differ from the previous
    base, so **every adjacent pair differs and the max homopolymer run is 1** (the
    strongest possible bound). A fixed virtual seed base defines the first base.
    ~1.58 bits/nt.
  - `encode(data: bytes, payload_bytes=16, index_trits=20, seed="A") -> (oligos, meta)`
  - `decode(oligos, meta) -> bytes`
  - `meta` carries `total_len`, `payload_bytes`, `index_trits`, `seed`. Passed
    out-of-band for M1 (a self-describing in-DNA header is a later enhancement; noted
    with a `ponytail:` marker in the code).

- `src/test_codec.py` - assert-based self-checks (no framework).

## Oligo format

Each oligo is one continuous rotating-coded base string over a trit stream:

```
[ index_trits (fixed width) | payload_trits (fixed width) ]
```

- Payload chunked into fixed `payload_bytes`; the last chunk is zero-padded and
  `decode` truncates to `meta.total_len`.
- `index`/`payload` integers are big-endian; `int <-> trits` is fixed-width base-3
  (smallest width `T` with `3^T >= 256^bytes`, checked with exact integer arithmetic).
- `decode` sorts oligos by decoded index, so reassembly is order-independent.

## Constraint story (kept honest)

- **Guaranteed:** max homopolymer run = 1, by construction.
- **Measured, not guaranteed:** GC content (reported as a distribution). A tunable
  hard GC window is the job of the screening codec (later), and the benchmark reports
  GC for both. M1 does not claim GC balancing.

## Tests (must pass)

- Exact round-trip on random payloads of length 0, 1, 15, 16, 17, 100, 1000.
- Every oligo has `max_homopolymer_run == 1`.
- Shuffled oligo list still decodes to the original (order independence).
- Index recovered correctly for every chunk.
- GC fraction reported and within (0, 1).

## Skipped for M1 (add when)

- ECC layer -> M2.
- Screening codec -> after M1 lands.
- Channel simulation and cited error rates -> M3.
- Self-describing in-DNA header -> when the pool must be parameter-free.
- Motif avoidance -> when the error model names motifs.
