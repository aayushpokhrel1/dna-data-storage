# M2 - Error correction (design)

Date: 2026-09-11. Status: done.

## Goal

A coding layer that recovers the original bytes after the errors a DNA channel
introduces (substitutions and, at M3, insertions/deletions), wrapping the M1
codec's `bytes <-> oligos` interface without changing its encoding.

## Decision (with the user, 2026-09-11)

Reed-Solomon + erasure now (with the rotating codec); fountain/LT codes added
later so the benchmark compares RS-vs-fountain as well as rotating-vs-screening
("both eventually, RS first").

## Model: detect and erase

- **Inner (detection):** each record carries a 4-byte CRC32 over `(global index +
  data)`. Any damaged oligo, substitution, garble, or a corrupted index, fails the
  CRC and is discarded. A corrupted index therefore cannot pollute another slot.
- **Outer (recovery):** systematic striped Reed-Solomon over GF(256) via `reedsolo`,
  one codeword per byte-column across the K data records of a block, adding
  `parity_records` parity records. A missing index (dropout) or CRC-failed oligo is
  an erasure; up to `parity_records` erasures per block are reconstructed.
- **Substitutions and indels both reduce to "bad oligo -> erasure."** An M3
  de-synchronized (indel) read fails its CRC and is recovered as an erasure.
- **Blocking:** data is split into `data_bytes`-sized records, grouped into blocks of
  `block_records` (with `block_records + parity_records <= 255`, the GF(256) limit),
  each block gets its own parity.

## Why not in-place inner correction

The M1 codec encodes a whole record as one base-3 big integer, so a single trit
error propagates through base-256 carries and can corrupt many bytes: byte-level
in-place correction is a poor fit. Detect-and-erase sidesteps this and keeps M1's
density. In-place correction (saving an erasure slot per substitution) needs
error-localized encoding (per-byte trits) and is a documented later enhancement.

## Gotchas recorded

- **Rotating-code error propagation:** an isolated substitution corrupts exactly two
  trits (the flipped position and the next, whose decode uses the corrupted base as
  its predecessor). Bounded, but it is why detection, not naive per-base handling, is
  the right inner layer. `decode` was made corruption-tolerant (a base that violates
  the rotating invariant emits a placeholder trit instead of raising).
- **Homopolymer runs are per-oligo:** each oligo is one physical molecule. Measuring
  `max_homopolymer_run` over a joined pool invents runs at oligo boundaries. The
  benchmark must measure per-oligo (the M1 guarantee is run = 1 per oligo).

## Interface

- `ecc.encode(data, data_bytes=8, parity_records=4, block_records=64, ...) -> (oligos, meta)`
- `ecc.decode(oligos, meta) -> bytes`  (raises `ecc.RecoveryError` past the budget)
- codec addition: `codec.decode_records(oligos, meta) -> {index: bytes}` returns only
  survivors so the ECC layer sees the erasures.

## Tests (`src/test_ecc.py`, all passing)

Clean round-trip; dropped oligos recovered; corrupted (payload/index) oligos
recovered; mixed drops + corruptions up to the budget recovered; exceeding the
budget raises `RecoveryError`; `decode_records` reports gaps.

## Skipped for M2 (add when)

- Fountain/LT codec -> after RS, for the benchmark comparison.
- In-place inner substitution correction -> if the recovery-vs-redundancy curve needs it.
- Channel simulation (real sub/ins/del rates, cited) -> M3.
