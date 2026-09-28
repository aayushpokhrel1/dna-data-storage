# Operations: running, verifying, regenerating

How to run this repo, which command produces each headline number, and the environment
traps that have cost time. Numbers here were re-verified on 2026-09-27.

Setup and the run commands live in the README's Reproducibility section, so they are
not repeated here. This file holds what the README does not: where each number comes
from, and the traps.

## What passing looks like

- `for t in src/test_*.py; do python "$t"; done` -> all nine exit 0, **53 checks**
  total (benchmark 12, ecc 8, channel 7, codec 6, screen_codec 6, fountain 5, pool 4,
  markercode 3, crossfamily 2).
- `python src/check_citations.py` -> 6 bib entries, 6 `\cite` keys, all defined.

## The operating point, and where each headline number comes from

Headline runs use `data_bytes=32` (~150-200 nt oligos), chosen because it is a
realistic synthesis length. Earlier numbers in this project's history were taken at
`data_bytes=16`, which flattered density and pushed the recovery cliff right; **any
number quoted without its operating point is not comparable.**

| Number | Value | Produced by |
| --- | --- | --- |
| Density, rotating+RS, 1x redundancy, coverage 8 | 0.634 bits/nt | `benchmark.compare_families` -> `results/family_comparison.json` |
| Density, screening+RS, 1x redundancy, coverage 8 | 0.780 bits/nt | same |
| Full-recovery threshold, all families except rotating+fountain | 1% per-base | same (rotating+fountain: 0.5%) |
| Density at low (4%) redundancy, rotating | 1.216 bits/nt | `src/scale_demo.py` (92,114 B, 2,999 oligos) |
| Density at low (4%) redundancy, screening | 1.498 bits/nt | `ecc.encode(..., base=screen_codec)` at the scale-demo settings (`data_bytes=32, parity_records=8, block_records=200`) |
| Coverage/parity frontier, rate 0.006 | full recovery at (coverage 8, 0.25x), (4, 0.5x), (2, 2x) | `benchmark.cost_grid` -> `results/cost_study.json` |
| RS vs fountain at 30% dropout, K=188 | RS full at 0.75x, fountain at 1.0x | `benchmark.dropout_study` -> same file |
| Indel confinement under deletions | ~0.84 of bases kept at 1%, ~0.94 at 0.5% (vs ~0.34 / ~0.45 unmarked) | `benchmark.indel_confinement` -> `indel.png` |

The screening density at low overhead is the one figure with no committed JSON behind
it; it is a one-line call, recorded in the table above so the paper's `~1.5 bits/nt`
claim is reproducible.

## Environment traps

- **Do not drive scripted edits to `paper/main.tex` through a shell heredoc.** The Bash
  tool collapses `\\` to `\`, which once turned `\num` into a literal newline plus `um`.
  Write the fix script to a file and run `python <file>`. (Also warned at the top of
  `main.tex` itself.)
- **Vault writes go through WSL, and note content must never sit on the command line**
  (the wrapper is re-quoted, so backticks and apostrophes get executed or break the
  command). Write the note to a file, then pass only literal paths.
- **Line endings.** The repo is on Windows; git warns LF -> CRLF on commit. Harmless.
- **Homopolymer runs are a per-oligo property.** Never measure `max_homopolymer_run`
  over a joined pool: it invents runs across oligo boundaries. `benchmark` already does
  this correctly and its docstring says so.
- **`reedsolo` works over GF(256)**, so an RS block cannot exceed 255 records. Parity is
  scaled to the block size, not to global K, in `benchmark._enc_kwargs`.

## Standing conventions

- **Cite, do not remember.** Every prior-art reference and every assumed error rate
  enters `paper/references.bib` only after being verified against the real source.
  All 6 entries were verified against the publisher or PubMed record on 2026-09-11.
- **No em dashes or en dashes** anywhere, per Aayush's global rule. Commas, colons,
  parentheses, or two sentences. In LaTeX this means `Reed-Solomon`, not `Reed--Solomon`.
- Claims stay computational and honest. This is a codec-and-benchmark paper, not a
  "DNA beats silicon" claim. Ours are modeled, the prior wet-lab numbers are measured,
  and the two are kept visibly separate.
- Non-trivial logic (encoder, decoder, ECC, channel) leaves a runnable check behind:
  round-trip and recovery assertions, not a demo.
- When referencing agentic tooling in any writeup, name **Claude Code** explicitly.
