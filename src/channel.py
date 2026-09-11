"""M3 synthesis/sequencing error channel for DNA-storage oligos.

A pure, deterministic function over a list of oligos: given a seed, the same
inputs always produce the same reads. Each oligo is a molecule; the channel
emits `coverage` reads of it (or none, if the molecule is dropped). Per read we
walk the oligo base by base and apply, independently at each position:

  - deletion  (p_del): the base is not emitted;
  - substitution (p_sub): a base drawn uniformly from the other three is emitted;
  - otherwise the base is emitted unchanged;
  - insertion (p_ins): after handling the base, one extra base drawn uniformly
    from all four is emitted.

Rates are parameters, not baked-in constants. The defaults are illustrative
only; no rate here is a literature value (PENDING CITATION, see the M3 spec).
"""
import os
import random
import sys

sys.path.insert(0, os.path.dirname(__file__))

BASES = "ACGT"


def _corrupt_read(oligo, rng, p_sub, p_ins, p_del):
    """Apply the per-base error model to one read of `oligo`."""
    out = []
    for base in oligo:
        if rng.random() < p_del:
            pass                                   # base deleted: emit nothing
        elif rng.random() < p_sub:
            alts = [b for b in BASES if b != base]
            out.append(rng.choice(alts))           # substitution: one of the other three
        else:
            out.append(base)                       # clean
        if rng.random() < p_ins:
            out.append(rng.choice(BASES))          # insertion: any of the four
    return "".join(out)


def corrupt(oligos, p_sub=0.0, p_ins=0.0, p_del=0.0, p_drop=0.0,
            coverage=1, seed=0):
    """Simulate reads of `oligos` under a seeded error channel.

    Returns a list of reads (str over "ACGT"). With probability `p_drop` an
    oligo is dropped entirely (emits no reads); otherwise it emits `coverage`
    independent reads. Deterministic for a given seed and inputs.
    """
    rng = random.Random(seed)
    reads = []
    for oligo in oligos:
        if rng.random() < p_drop:
            continue                               # molecule lost
        for _ in range(coverage):
            reads.append(_corrupt_read(oligo, rng, p_sub, p_ins, p_del))
    return reads


if __name__ == "__main__":
    pool = ["ACGTACGTAC", "GTACGTACGT"]
    print("clean  :", corrupt(pool, seed=0))
    print("noisy  :", corrupt(pool, p_sub=0.1, p_ins=0.05, p_del=0.05, seed=7))
