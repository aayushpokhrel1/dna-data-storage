# M3 - Channel and benchmark (design + delegation spec)

Date: 2026-09-11. Status: implementing (tests in-house, implementation delegated).

This doc is also the self-contained spec for the delegated worker. It has no
conversation context: build `src/channel.py` and `src/benchmark.py` to satisfy the
tests in `src/test_channel.py` and `src/test_benchmark.py` and the interfaces below.
Do not modify `src/codec.py`, `src/constraints.py`, `src/ecc.py`, or any test file.

## Goal

Simulate synthesis/sequencing errors, run data through encode -> corrupt -> decode,
and measure recovery vs error rate, density, and overhead. The recovery-vs-rate
sweep is the deliverable; error rates are parameters, not baked-in constants.

## Existing interfaces to build on (do not change)

- `ecc.encode(data, data_bytes=8, parity_records=4, block_records=64, index_trits=20, seed="A") -> (oligos, meta)`
- `ecc.decode(oligos, meta) -> bytes` (raises `ecc.RecoveryError` past the parity budget)
- `constraints.gc_fraction(seq) -> float`, `constraints.max_homopolymer_run(seq) -> int`
- oligos are `str` over "ACGT"; a "read" is also such a string.

## `src/channel.py`

Deterministic given a seed. A pure function over a list of oligos.

- `corrupt(oligos, p_sub=0.0, p_ins=0.0, p_del=0.0, p_drop=0.0, coverage=1, seed=0) -> list[str]`
  - For each oligo, emit `coverage` reads (unless dropped). With probability `p_drop`
    the oligo is dropped entirely (emits no reads).
  - For each read, walk the oligo base by base:
    - with prob `p_del`, delete this base (emit nothing);
    - else with prob `p_sub`, emit a base chosen uniformly from the other three;
    - else emit the base unchanged;
    - after handling the base, with prob `p_ins` emit one extra base chosen uniformly
      from all four (an insertion).
  - Order of the returned reads must not encode the original index (shuffle within a
    fixed seed is fine, or keep order; the decoder is order-independent either way).
  - Same seed + same inputs => identical output (use a local `random.Random(seed)`).
- Illustrative defaults only; PENDING CITATION. Do not present any rate as a
  literature value. A realistic operating point is cited at the paper stage.

## `src/benchmark.py`

- `code_rate(data_len, oligos, meta) -> dict` with at least:
  - `bits_per_nt = data_len * 8 / total_nt` where `total_nt = sum(len(o) for o in oligos)`
  - `code_rate = data_len / (data_len + parity_bytes_equivalent)` (report the parity
    fraction honestly; document the formula you use in a comment)
  - `n_oligos`, `total_nt`
- `constraint_stats(oligos) -> dict`: per-oligo GC (min/mean/max) and per-oligo
  `max_homopolymer_run` (max over oligos). NEVER measure a run over joined oligos.
- `recovery_at_rate(data, rate, trials, seed, ratio=(1,1,1), **ecc_kwargs) -> float`:
  split `rate` across (sub, ins, del) by `ratio`, run `trials` independent
  encode->corrupt->decode passes (vary the channel seed per trial), return the
  fraction that decode exactly to `data`. A `RecoveryError` or a mismatch counts as a
  failure (must not raise out of this function).
- `sweep(data, rates, trials, seed, **kwargs) -> dict`: run `recovery_at_rate` over
  `rates`, return `{"rates": [...], "recovery": [...], "config": {...}}`.
- `main(smoke=False)`: run a small sweep (smoke: tiny/fast), write
  `results/benchmark.json` with the sweep result, `code_rate`, and `constraint_stats`,
  and print a short summary. Guard with `if __name__ == "__main__":` and an
  `argparse` `--smoke` flag.

## Metrics (definitions fixed here, in-house)

- **Recovery** at a rate = fraction of trials that decode exactly to the input.
- **Max recoverable rate** = highest swept rate with recovery == 1.0.
- **Density** = payload bits / total encoded nt (bits/nt).
- **Overhead / code rate** = payload bytes / encoded record bytes (parity + CRC +
  index all count as overhead; document the formula).

## Self-checks (the tests encode these; all must pass)

Channel: zero rates => identity; `p_sub=1` preserves length and changes every base;
`p_del=1` empties reads; `p_ins=1` lengthens; `p_drop=1` drops all; determinism by seed.
Benchmark: recovery == 1.0 at rate 0; recovery == 1.0 for a low rate within the parity
budget; recovery < 1.0 (no crash) at a high rate; `sweep` returns aligned rate/recovery
lists; `code_rate`/`constraint_stats` return the documented keys with per-oligo runs.

## Landed beyond the first cut

- **Sequencing coverage** turned out to be necessary, not optional: with a single
  read, ~81 nt oligos at a 1% per-base rate lose ~55% of oligos, so the outer code
  cliffs at <0.5% and the curve is uncomparable to the coverage-based literature.
  `ecc.decode` now keeps the first CRC-valid read per index, so `coverage` reads let
  one clean read fill a slot. With `coverage=10` the curve is realistic: full
  recovery through ~1%, ~0.85 at 2%, failing by 5% (`results/benchmark.json`).
- **Cited error context (verified 2026-09-11, `paper/references.bib`):** the swept
  range and the sub:ins:del ratio (5:1:1) are framed on Xu et al. 2021 (synthesis
  ~0.7%: sub 0.5 / ins 0.1 / del 0.1; Illumina <1%; nanopore ~10%). Codec prior art:
  Church 2012, Goldman 2013, Grass 2015, Erlich 2017, Organick 2018. `src/check_citations.py`
  validates bib integrity and (later) \cite coverage.

## Skipped for M3 (add when)

- Polished publication figures -> M4.
- Majority-vote consensus across coverage reads (vs the current first-CRC-valid-read)
  -> if the recovery curve needs squeezing further.
- Density sweep over `payload_bytes`/`index_trits` (M1 note) -> M4 density figure;
  the current 0.70 bits/nt reflects the small `data_bytes=8` default, index/CRC-dominated.
