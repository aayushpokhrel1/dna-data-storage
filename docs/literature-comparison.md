# Literature context

Where this codec sits relative to the canonical DNA data-storage works. Every entry
is verified against its source (see `paper/references.bib`, verified 2026-09-13).

**Read this honestly.** The prior works are **wet-lab demonstrations**: they
synthesized real DNA, sequenced it, and recovered real files. This work is a
**computational framework**: encode, a simulated synthesis/sequencing channel, decode,
measured in software. So the numbers are not directly comparable, ours are modeled,
theirs are measured. The contribution here is not a density record but a unified,
reproducible comparison of codec and error-correction families under one channel, plus
random access and an integrated indel layer. Density is quoted only where a specific
bits/nt figure was verified against the source; qualitative cells avoid inventing one.

| Work | Milestone (measured) | Error correction | Indel handling | Random access | Net density (bits/nt) |
| --- | --- | --- | --- | --- | --- |
| Church 2012 [`church2012`] | ~5.3 Mbit encoded in oligos | none (small bit-error rate accepted) | none | no | ~1 gross, no ECC (qualitative) |
| Goldman 2013 [`goldman2013`] | 739 KB, recovered at 100% | fourfold overlapping redundancy | via redundancy | no | low, 4x redundancy (qualitative) |
| Grass 2015 [`grass2015`] | file in silica, error-free | Reed-Solomon (inner + outer) | as erasures | no | **1.14** (verified) |
| Erlich 2017 [`erlich2017`] | 2.1 MB, near Shannon capacity | fountain (LT) + RS droplet screen | dropouts as erasures | no | **1.57** (coding potential 1.98) (verified) |
| Organick 2018 [`organick2018`] | 200 MB, 13M oligos, error-free | RS + read consensus | as erasures | yes (PCR primers) | competitive (qualitative) |
| **This work (modeled)** | framework + benchmark, in simulation | **RS and fountain (both, compared)** | **marker resync + inner RS (integrated)** | yes (file-id barcode) | 0.780 at 1x redundancy; 1.498 at 4% (screening), 1.216 (rotating) |

## What the comparison says

- **Density.** Our screening codec reaches 1.498 bits/nt at low (4%) redundancy and
  realistic ~150-200 nt oligos, between Grass's Reed-Solomon 1.14 and Erlich's fountain
  1.57, which is what one expects from using the same coding families. The rotating codec
  trades density (1.216 bits/nt at the same settings) for a homopolymer-run-of-1
  guarantee. Both are measured at the scale-demo settings (`data_bytes=32`,
  `parity_records=8`, `block_records=200`) on a 92,114-byte payload in 2,999 oligos; the
  rotating figure is the committed `src/scale_demo.py` run and the screening figure is the
  same call with `base=screen_codec` (see `docs/OPERATIONS.md`). We do not claim a density
  record; the two verified anchors (Grass, Erlich) frame the range.
- **Coverage.** Erlich recovered from ~1.3% oligo dropout and Organick from ~10% nanopore
  error, both relying on sequencing coverage. Our cost study reproduces that dependence:
  coverage is the dominant recovery lever (below coverage 2 nothing recovers at rate
  0.006, whatever the parity), and coverage and parity are substitutable budgets above
  that.
- **What's new here.** No prior work provides an open, apples-to-apples comparison of
  {rotating, screening} x {RS, fountain} under one channel, and the honest finding that
  RS (MDS) matches or beats fountain at fixed redundancy (fountain's edge is
  ratelessness). We also integrate a marker-resync + inner-RS indel layer and show it
  only works with an error-localized base codec, a point implicit in the prior RS-erasure
  designs but not, to our knowledge, measured this way.
- **Random access.** Organick demonstrated wet-lab random access with PCR primers; we
  model the codec-level analogue (a file-id barcode) and show per-file retrieval is O(1)
  in the archive size.

## Verification notes

- Grass 2015 net density 1.14 bits/nt and Erlich 2017 net density 1.57 bits/nt (coding
  potential 1.98) are stated in Erlich & Zielinski 2017 and confirmed against the paper.
- Church (~5.3 Mbit), Goldman (739 KB, 100% recovery, fourfold redundancy), and Organick
  (200 MB, 13M oligos, error-free, random access) figures are confirmed against their
  sources; their exact net bits/nt were not independently re-derived here, so the density
  cell is left qualitative rather than quote an unverified number.
