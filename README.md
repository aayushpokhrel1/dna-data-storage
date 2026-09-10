# DNA Data Storage: A Constraint-Aware Codec with Error Correction

> **Status: M1 done.** The constraint-aware codec (rotating code) encodes bytes to
> DNA oligos and back, exact round-trip, homopolymer runs bounded to 1 by
> construction, indexed oligos for an addressable pool. Code in `src/codec.py` and
> `src/constraints.py`, checks in `src/test_codec.py`. Later sections are still
> marked **planned** so it is always clear what exists and what does not.

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
- [ ] **M2 - Error correction.** Add the coding layer; recover from substitutions
      and from insertions and deletions.
- [ ] **M3 - Channel and benchmark.** Simulate synthesis and sequencing errors;
      measure density, overhead, and recovery versus error rate.
- [ ] **M4 - Figures.** Recovery curves, density and overhead, constraint
      distributions.
- [ ] **M5 - Paper draft.** Write the manuscript in `paper/`.
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
