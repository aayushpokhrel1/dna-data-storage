# DNA Data Storage: A Constraint-Aware Codec with Error Correction

> **Status: M4 done.** Constraint-aware rotating codec (M1), Reed-Solomon error
> correction (M2), a seeded channel + recovery benchmark (M3), and publication
> figures (M4). The codec recovers fully through ~1% per-base error (coverage 10),
> and density rises from 0.70 toward the rotating ceiling (log2 3) as the payload
> grows. Figures in `results/figures/`, error rates and prior art cited in
> `paper/references.bib` (verified 2026-09-11). A substance program (S1-S6) is now
> underway to turn the tool into a research contribution, an open framework comparing
> codec/ECC families across density, robustness, and cost: the screening codec (S1,
> `src/screen_codec.py`, ~2 bits/nt with a guaranteed GC window) and the fountain
> codec (S2, `src/fountain.py`), and the cross-family benchmark + density-robustness
> Pareto frontier (S3), the coverage/parity cost study (S4), random access (S5), and a marker resync inner
> code for indels (S6) are done, the full substance program. screening+RS dominates the
> density-robustness frontier, RS matches or beats fountain at fixed redundancy
> (fountain's edge is ratelessness), per-file retrieval is O(1) in the archive size,
> and marker resync + inner RS (integrated into the pipeline) recovers deletions the
> detect-and-erase baseline loses. Headline runs use a realistic ~150-200 nt oligo
> length (screening+RS ~0.78 bits/nt at 1x redundancy). Code in `src/` with 52 passing
> assertion checks across nine `test_*.py` files. The paper draft (M5) is next.

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
      robustness (`results/figures/pareto.png`). Finding: screening lifts density
      (~0.53 -> 0.64 bits/nt at 1x overhead) and screening+RS dominates the frontier
      here; fountain trails RS at this small block size (its edge is large,
      dropout-heavy pools). Checks in `src/test_crossfamily.py` and
      `src/test_benchmark.py`.
- [x] **S4 - Coverage/parity cost study.** `benchmark.cost_grid` sweeps recovery
      over the (coverage x redundancy) grid at a fixed error rate ->
      `results/figures/cost_grid.png` with an iso-recovery frontier: coverage and
      parity are substitutable budgets (full recovery at coverage 16 / 0.25x
      overhead == coverage 8 / 1x == coverage 4 / 2x). `benchmark.dropout_study`
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

- [ ] **M5 - Paper draft.** Write the manuscript in `paper/`, covering the framework,
      the codec/ECC families, the comparison, random access, and the cost study.
- [ ] **M6 - Preprint.**

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
├── results/                  # metrics JSON and figures
└── paper/                    # references.bib (+ LaTeX draft, planned)
```

## License

- **Code** under the MIT License (see [`LICENSE`](LICENSE)).
- **Paper text and figures** under CC-BY 4.0 (see [`LICENSE-paper`](LICENSE-paper)).
