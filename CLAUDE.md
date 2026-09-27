# Working in this repo

DNA data storage research (Python).

## Where things are written down

| File | What it holds |
| --- | --- |
| `README.md` | What it is and how to run it |
| `docs/literature-comparison.md` | Comparison against the literature |
| `docs/superpowers/` | Per-feature specs and plans |

Python; dependencies in `requirements.txt`. For research, a **result and the method that produced it** belong in the repo docs, while the *finding* and what it changes belong in the vault.

## The handover

`HANDOVER.md` holds **current state only**: where things stand, what is half-done, what is next.

**Before adding a line to it, ask: will this still be true in a month?** If yes it belongs in one
of the permanent docs above, or in the vault. This rule exists because the same file on a sibling
project reached 1,365 lines by accumulating everything and began contradicting itself, which
caused real wrong work: a session acted on a superseded line and re-proposed something the same
file recorded as already tried and rejected.

`HANDOVER.md` is **gitignored**, so it is local to this machine and has no git history. Nothing
durable can survive there, and it cannot be recovered if overwritten - copy it to a scratchpad
before any rewrite.

## When Aayush says "update"

"Update the docs", "update everything", or just "update" means **all of it, in this turn**:

1. **The knowledge vault** - `C:\Users\aayus\Documents\Knowledge-Vault\Projects\dna-data-storage\`
   (the path is per-machine). Add what this session learned that is worth keeping: a decision and
   its WHY, a non-obvious gotcha or fix, a research finding, a cross-project learning. This is the
   part that gets forgotten, and it is the part that compounds.
2. **Every doc in this repo**, not only the one already open.
3. **The handover**, but only its current-state numbers.

**Vault writes go through WSL, and note content must never appear on the command line.** The Bash
tool re-quotes the wrapper, so backticks and apostrophes inside a note get executed or break the
command. Write the note to a file first, then pass only literal paths:

```
wsl -d Ubuntu -- bash -lc 'cat /mnt/c/<tmp>/note.md >> /mnt/c/Users/aayus/Documents/Knowledge-Vault/Projects/dna-data-storage/index.md'
```

**Updating docs means making them TRUE, not just appending what shipped.** Correct or strike a
stale claim where it sits rather than adding a newer entry underneath it, because the next reader
may hit the old one first. Cross-check every number (versions, counts, commit) against reality
instead of trusting what the file says.

## Where knowledge goes

A lesson has exactly one home, chosen by how far it reaches:

- **A rule about specific code goes in a comment AT that code.** The most reliable form there is:
  you cannot edit the function without reading the warning above it. A note filed elsewhere is the
  least reliable, because nobody goes back to read it.
- **A lesson that generalises goes in the vault**, phrased so it is useful on a different
  project, with this one as the example.
- **A dated narrative of what you did today goes in `git log`.** It is already there, in detail.

Keep any one lesson in a single place. Two copies drift, and the drift is what causes wrong work
later.
