# DNA Data Storage: A Constraint-Aware Codec with Error Correction

> **Status: M4 done.** Constraint-aware rotating codec (M1), Reed-Solomon error
> correction (M2), a seeded channel + recovery benchmark (M3), and publication
> figures (M4). The codec recovers fully through ~1% per-base error (coverage 10),
> and density rises from 0.70 toward the rotating ceiling (log2 3) as the payload
> grows. Figures in `results/figures/`, error rates and prior art cited in
> `paper/references.bib` (verified 2026-09-11). A substance program (S1-S6) is now
> underway to turn the tool into a research contribution, an open framework comparing
> codec/ECC families across density, robustness, and cost: the screening codec (S1,
> `src/screen_codec.py`, ~2 bits/nt with a guaranteed GC window) is done. Code in
> `src/` with 33 passing assertion checks across six `test_*.py` files. Remaining
> sections are marked **planned**.

Working title: *A Constraint-Aware Codec for DNA Data Storage: Encoding, Error
Correction, and a Recovery Benchmark* (not final).

This project builds and benchmarks a codec that stores digital data in DNA:
mapping bits to ACGT sequences under real biological constraints, adding error
correction, simulating the synthesis-and-sequencing error channel, and measuring
how much data survives. The working tool and its measured performance are the
contribution.

## Scope

This paper is about **DNA data storage**: encoding digital bits in synthesized DNA
sequences that are later read by sequencing. It is a sequence-level, information-
theoretic problem and is **independent of the charge-transport physics** in the
first two projects. It does not reuse the paper-1 solver.

- Paper 1: `aayushpokhrel1/dna-quantum-surrogates` (DNA hole-transport simulation).
- Paper 2: `aayushpokhrel1/dna-vs-silicon-memory` (molecular-electronic memory).
- Paper 3 (this repo): DNA data storage codec. Separate field, separate repo.

## Contribution (planned)

A runnable, tested codec plus a benchmark, kept honest and measurable:

1. **Constraint-aware encoding.** Map arbitrary bytes to DNA under the constraints
   real synthesis and sequencing impose (balanced GC content, no long homopolymer
   runs, avoid problematic motifs).
2. **Error correction.** Add a coding layer that recovers the original data after
   the errors DNA channels introduce, including the hard case of insertions and
   deletions, not just substitutions.
3. **Channel simulation and benchmark.** Model synthesis and sequencing errors at
   configurable rates, run data through encode, corrupt, decode, and report the
   metrics below.

## Metrics (planned)

- **Information density:** bits stored per nucleotide (theoretical maximum is 2;
  constraints and error-correction overhead reduce it).
- **Coding overhead:** redundancy added by the error-correction layer.
- **Recovery:** fraction of data recovered as a function of the channel error rate,
  and the maximum error rate that still recovers everything.
- **Constraint satisfaction:** measured GC content and homopolymer-run distribution
  of the encoded sequences.
- **Context:** the theoretical density of DNA storage next to silicon and flash, as
  illustrative background, not a headline claim.

## Open questions (to decide together before building)

- **Which error-correction scheme:** Reed-Solomon (handles substitutions well),
  fountain and LT codes (used by high-density DNA-storage work), or LDPC. Insertions
  and deletions need extra handling (synchronization or marker codes), so the indel
  strategy is a key decision.
- **Error model:** what substitution, insertion, and deletion rates to assume for
  synthesis and sequencing, and from which source.
- **Random access / addressing:** whether the codec needs to retrieve a single file
  from a pool, or only encode and decode a whole payload.
- **How far the silicon comparison goes:** density only, or also cost, latency, and
  retention.

## Roadmap (planned)

- [x] **M1 - Encoding.** Constraint-aware rotating-code codec in `src/codec.py`
      (bytes <-> DNA oligos, ECC-agnostic) with shared constraint measurement in
      `src/constraints.py`. Homopolymer runs are bounded to 1 by construction and the
      codec round-trips exactly with no channel; GC is measured and reported (a hard
      GC window is the screening codec's job, later). Oligos are indexed so decode
      reassembles from a shuffled or partial pool. Checks in `src/test_codec.py`.
      Design: `docs/superpowers/specs/2026-09-10-m1-encoding-design.md`.
- [x] **M2 - Error correction.** Reed-Solomon layer in `src/ecc.py` wrapping the M1
      codec, detect-and-erase: a per-oligo CRC (over index + data) flags any damaged
      oligo, and a systematic striped RS outer code (over GF(256), via `reedsolo`)
      reconstructs up to `parity_records` erased or corrupted oligos per block.
      Substitutions and indels both reduce to "bad oligo -> erasure". Checks in
      `src/test_ecc.py`. Design:
      `docs/superpowers/specs/2026-09-11-m2-ecc-design.md`.
- [x] **M2b - Screening codec.** Higher-density second codec family in
      `src/screen_codec.py` (DNA Fountain-style, Erlich 2017): bits pack 2/nt and a
      per-chunk seed is screened by rejection so the oligo satisfies BOTH a GC window
      and the homopolymer bound by construction (the GC guarantee the rotating codec
      lacks). ~2 bits/nt at large payloads. Same indexed, ECC-agnostic interface as
      the rotating codec. Checks in `src/test_screen_codec.py`.
- [x] **M3 - Channel and benchmark.** Seeded synthesis/sequencing channel
      (`src/channel.py`: per-base substitution/insertion/deletion, whole-oligo
      dropout, sequencing coverage) and a recovery benchmark (`src/benchmark.py`:
      recovery vs per-base rate, density in bits/nt, code rate, per-oligo GC and
      homopolymer stats -> `results/benchmark.json`). Error range and sub:ins:del
      ratio framed on cited literature (`paper/references.bib`, verified against the
      sources); `src/check_citations.py` guards citation integrity. Checks in
      `src/test_channel.py` and `src/test_benchmark.py`. Design:
      `docs/superpowers/specs/2026-09-11-m3-channel-benchmark-design.md`.
- [x] **M4 - Figures.** `src/figures.py` -> `results/figures/`: recovery vs error
      rate (one line per sequencing coverage, showing coverage as the recovery
      lever), information density vs payload size (against the log2 3 rotating
      ceiling, from `benchmark.density_sweep` with the data-record count held fixed
      so the parity fraction stays constant), and per-oligo GC + homopolymer-run
      distributions (all runs = 1, the M1 guarantee). Plotting delegated to a
      cheap-model worker and reviewed; the density-sweep computation is in-house
      (checked in `src/test_benchmark.py`).
### Substance program (approved 2026-09-12): turn the tool into a contribution

The core new claim is an open, reproducible framework for comparing constraint-aware
DNA-storage codecs across the density / robustness / cost trade-off, with random
access. Kept computational and honest, no "DNA beats silicon".

- [x] **S1 - Screening codec.** See M2b above (`src/screen_codec.py`).
- [ ] **S2 - Fountain codec.** LT/fountain outer code (Erlich 2017) as a second ECC
      family, so the benchmark compares RS-vs-fountain.
- [ ] **S3 - Cross-family benchmark + Pareto.** One harness over {rotating, screening}
      x {RS, fountain}; report the density-robustness Pareto frontier.
- [ ] **S4 - Coverage/parity cost study.** Pareto-optimal (coverage x redundancy)
      budget for a target recovery at a given error rate.
- [ ] **S5 - Random-access demo.** Retrieve one file from a mixed pool of many;
      measure reads needed. Realizes the addressable-oligo design.
- [ ] **S6 - Indel-correcting inner code.** Marker/watermark code that corrects
      insertions/deletions in place instead of discarding oligos (the stretch novelty).

- [ ] **M5 - Paper draft.** Write the manuscript in `paper/`, covering the framework,
      the codec/ECC families, the comparison, random access, and the cost study.
- [ ] **M6 - Preprint.**

## Repository layout

```
dna-data-storage/
├── README.md          # this file
├── LICENSE            # MIT (code)
├── LICENSE-paper      # CC-BY 4.0 (paper text and figures)
├── requirements.txt   # Python dependencies
├── src/               # codec, channel, and benchmark code
├── data/              # small or synthetic inputs (large data is gitignored)
├── results/           # figures and metrics
└── paper/             # LaTeX draft
```

## License

- **Code** under the MIT License (see [`LICENSE`](LICENSE)).
- **Paper text and figures** under CC-BY 4.0 (see [`LICENSE-paper`](LICENSE-paper)).
