# DNA Data Storage: A Constraint-Aware Codec with Error Correction

> **Status: M5 done, manuscript drafted. Next is M6 (author review, then arXiv).**
> The full pipeline is built and measured: two constraint-aware base codecs (rotating,
> screening), two error-correction families (Reed-Solomon detect-and-erase, LT
> fountain), a seeded synthesis/sequencing channel, an integrated marker-resync plus
> inner-RS indel layer, codec-level random access, and a cross-family benchmark, with
> the manuscript in [`paper/main.tex`](paper/main.tex) (self-contained, compiles on
> Overleaf/arXiv).
>
> Headline numbers, all at the realistic operating point (`data_bytes=32`, ~150-200 nt
> oligos): screening lifts density from 0.634 to 0.780 bits/nt at 1x redundancy and
> coverage 8, and reaches 1.498 bits/nt at 4% redundancy (rotating: 1.216, measured on
> a 92 kB real-file run). Three of the four families hold full recovery to 1% per-base
> error; only rotating+fountain trails, at 0.5%. Coverage and redundancy are
> substitutable budgets; RS matches or beats fountain at fixed redundancy, so
> fountain's edge is ratelessness, not erasure efficiency; per-file retrieval is O(1)
> in archive size; and the indel layer only works with an error-localized base codec.
> Code in `src/` with 53 passing assertion checks across nine `test_*.py` files.
> Prior art and error rates cited in `paper/references.bib` (verified 2026-09-11).
> Which command produces each number: [`docs/OPERATIONS.md`](docs/OPERATIONS.md).

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

## Contribution

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

## Metrics

- **Information density:** bits stored per nucleotide (theoretical maximum is 2;
  constraints and error-correction overhead reduce it).
- **Coding overhead:** redundancy added by the error-correction layer.
- **Recovery:** fraction of data recovered as a function of the channel error rate,
  and the maximum error rate that still recovers everything.
- **Constraint satisfaction:** measured GC content and homopolymer-run distribution
  of the encoded sequences.
- **Context:** the theoretical density of DNA storage next to silicon and flash, as
  illustrative background, not a headline claim.

Positioning against the canonical DNA-storage works (Church, Goldman, Grass, Erlich,
Organick), with verified numbers and a clear modeled-vs-measured split, is in
[`docs/literature-comparison.md`](docs/literature-comparison.md).

## Settled design decisions

These were open questions early on and are now answered and built. They are recorded
here so they are not re-opened by accident.

- **Error correction: both families, compared.** Reed-Solomon (detect-and-erase, a
  per-oligo CRC turns any damage into an erasure) and an LT fountain code, benchmarked
  against each other rather than one chosen. LDPC was not pursued.
- **Indels: marker resync plus an inner RS code**, integrated into the pipeline and
  gated by `inner_nsym` / `marker_period`. This confines an indel's damage and corrects
  the residual local byte errors. Full in-place correction (Davey-MacKay watermark plus
  soft decoding) is deliberately out of scope and named as future work.
- **Error model: taken from the literature, not assumed.** Synthesis ~0.7% (sub 0.5 /
  ins 0.1 / del 0.1), Illumina <1%, nanopore ~10%, from Xu 2021 (`xu2021`), verified
  against the source. The benchmark uses a 5:1:1 sub:ins:del ratio.
- **Random access: yes, at the codec level.** A file-id barcode per oligo
  (`src/pool.py`); retrieval filters the pool and decodes one file's subset. Physical
  primer-based retrieval is out of scope; the barcode is the codec-level address.
- **Silicon comparison: illustrative context only.** Density background, never a
  headline claim. No "DNA beats silicon" framing anywhere.
- **Scope and licensing.** DNA data storage at the sequence level, independent of the
  charge-transport physics in papers 1 and 2. MIT for code, CC-BY 4.0 for the paper.

Still open, framed as future work in the draft rather than blockers: a rateless or
incremental experiment to show fountain's true advantage; wiring the marker layer into
the fountain path (currently RS-only); full Davey-MacKay indel correction; validation
against measured per-position error profiles; and majority-vote consensus decoding
instead of the current first-CRC-valid-read.

## Roadmap

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
- [x] **S2 - Fountain codec.** LT/fountain outer code in `src/fountain.py` (Erlich
      2017): data split into K segments, N>K droplets each XOR a Robust-Soliton subset
      seeded by the droplet number (= codec index), per-droplet CRC discards corrupt
      droplets, peeling decoder reconstructs. Dropouts and corruptions are both just
      missing droplets. Same `encode/decode` interface as `ecc.py`. Checks in
      `src/test_fountain.py`.
- [x] **S3 - Cross-family benchmark + Pareto.** ECC layers made codec-pluggable
      (`ecc`/`fountain` take a `base` codec); `benchmark.compare_families` sweeps all
      {rotating, screening} x {RS, fountain} at matched redundancy and coverage ->
      `results/family_comparison.json`, and `figures.fig_pareto` plots density vs
      robustness (`results/figures/pareto.png`). Finding at the current operating point
      (1x overhead, coverage 8): screening lifts density 0.634 -> 0.780 bits/nt, and
      both screening families share the frontier point (equal density, both holding full
      recovery to 1% per-base, as does rotating+RS at the lower density). Only
      rotating+fountain trails, at 0.5%. So the base codec, not the error-correction
      family, separates the four here. Checks in `src/test_crossfamily.py` and
      `src/test_benchmark.py`.
- [x] **S4 - Coverage/parity cost study.** `benchmark.cost_grid` sweeps recovery
      over the (coverage x redundancy) grid at a fixed error rate ->
      `results/figures/cost_grid.png` with an iso-recovery frontier: coverage and
      parity are substitutable budgets (at rate 0.006, full recovery at coverage 8 /
      0.25x overhead, coverage 4 / 0.5x, and coverage 2 / 2x). `benchmark.dropout_study`
      compares RS vs fountain at 30% whole-oligo dropout and large K
      (`results/figures/dropout.png`). Honest finding: RS (MDS-optimal for erasures)
      matches or beats fountain at fixed redundancy (full recovery at 0.75x vs 1x);
      fountain's real value is ratelessness (droplets on demand), which a fixed-N
      benchmark does not reward. Data in `results/cost_study.json`. Checks in
      `src/test_benchmark.py`.
- [x] **S5 - Random-access demo.** `src/pool.py` stores many files in one shuffled,
      barcoded oligo pool; retrieving a file filters the pool by its file-id barcode
      and decodes only that subset, reusing any base codec + ECC layer. Barcode
      errors just drop an oligo from its file (the ECC recovers it). Per-file
      retrieval reads one file's worth of oligos regardless of archive size (O(1) vs
      the O(M) whole-pool decode): `results/figures/random_access.png`. Checks in
      `src/test_pool.py` (roundtrip, subset-only access, cross-file isolation,
      recovery under channel).
- [x] **S6 - Indel handling: marker resync + inner RS, integrated.**
      `src/markercode.py` inserts a known 4-mer marker every `period` bases and
      re-anchors on decode, confining an indel to one run. This is now wired into the
      pipeline: `ecc.encode(..., inner_nsym=, marker_period=)` adds a per-oligo inner
      Reed-Solomon code so the residual local byte errors a resynced indel leaves are
      corrected, and the oligo is recovered instead of erased. End-to-end
      (`benchmark.indel_pipeline` -> `results/figures/indel.png`): with the screening
      base codec the integrated codec recovers deletions up to ~1% that the
      detect-and-erase baseline loses. Key finding: this requires an error-LOCALIZED
      base codec (screening's 2-bit packing) -- the rotating codec's big-integer
      records spread one base error across the whole record, so markers alone cannot
      help it (the M2 error-propagation point). Full in-place indel correction
      (Davey-MacKay watermark + soft decoding) remains out of scope. Checks in
      `src/test_markercode.py` and `src/test_ecc.py`.

- [x] **M5 - Paper draft.** Self-contained manuscript in [`paper/main.tex`](paper/main.tex)
      (compiles on Overleaf/arXiv, `pdflatex` + `bibtex`): the framework claim, the two
      base codecs and two ECC families, the channel and metrics, and the results with the
      eight real figures and the honest findings (screening+RS on the frontier; RS beats
      fountain at fixed redundancy; coverage/parity substitution; O(1) random access;
      indel correction needs error-localized coding). Every DNA-side number is from the
      tested pipeline; every prior-art number is verified (`paper/references.bib`,
      `docs/literature-comparison.md`).
- [ ] **M6 - Preprint.** Author review, then post to arXiv.

## Reproducibility

Pure Python, standard library plus `numpy`, `matplotlib`, and `reedsolo`.

```bash
python -m venv .venv
# activate, then:
pip install -r requirements.txt
```

Run the tests (assert-based self-checks, no framework):

```bash
for t in src/test_*.py; do python "$t"; done
```

Run the benchmark (recovery sweep, cross-family comparison, cost study) and write
the JSON results:

```bash
python src/benchmark.py            # writes results/*.json
python src/benchmark.py --smoke    # fast subset
```

Regenerate every figure from the measured numbers:

```bash
python src/figures.py              # writes results/figures/*.png
```

Run the end-to-end scale demo (stores this codec's own source, ~90 KB / ~3000
oligos, in DNA and recovers it exactly through a realistic channel):

```bash
python src/scale_demo.py
```

Check citation integrity (bib fields present, every `\cite` defined once a draft
exists):

```bash
python src/check_citations.py
```

## Repository layout

```
dna-data-storage/
├── README.md
├── LICENSE / LICENSE-paper   # MIT (code), CC-BY 4.0 (paper)
├── requirements.txt
├── src/
│   ├── constraints.py        # GC / homopolymer / motif measurement
│   ├── codec.py              # rotating base codec (homopolymer run = 1)
│   ├── screen_codec.py       # screening base codec (2 bits/nt, GC window)
│   ├── ecc.py                # Reed-Solomon detect-and-erase ECC layer
│   ├── fountain.py           # LT/fountain ECC layer
│   ├── markercode.py         # marker resync inner code (indels)
│   ├── channel.py            # synthesis/sequencing error channel
│   ├── pool.py               # random access over a mixed file pool
│   ├── scale_demo.py         # end-to-end ~90 KB real-file storage + recovery
│   ├── benchmark.py          # metrics, sweeps, cross-family + cost studies
│   ├── figures.py            # publication figures
│   ├── check_citations.py    # citation integrity
│   └── test_*.py             # per-module self-checks
├── data/                     # small or synthetic inputs (large data gitignored)
├── docs/
│   ├── OPERATIONS.md         # run, verify, which command makes which number, traps
│   ├── literature-comparison.md
│   └── superpowers/specs/    # per-feature design records
├── results/                  # metrics JSON and figures
└── paper/
    ├── main.tex              # manuscript draft
    └── references.bib        # verified citations
```

## License

- **Code** under the MIT License (see [`LICENSE`](LICENSE)).
- **Paper text and figures** under CC-BY 4.0 (see [`LICENSE-paper`](LICENSE-paper)).
