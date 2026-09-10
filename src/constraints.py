"""Biological-constraint measurement for DNA-storage sequences.

Pure predicates and measurement, shared by the codecs and the benchmark. No
encoding logic lives here.
"""
from itertools import groupby


def gc_fraction(seq):
    """Fraction of bases that are G or C. Empty sequence -> 0.0."""
    if not seq:
        return 0.0
    return sum(b in "GC" for b in seq) / len(seq)


def max_homopolymer_run(seq):
    """Length of the longest run of one identical base. Empty -> 0."""
    return max((len(list(g)) for _, g in groupby(seq)), default=0)


def has_forbidden_motif(seq, motifs=()):
    """True if any forbidden motif occurs in seq. Stub until the error model
    (M3) names real motifs; the interface is here so callers do not change."""
    return any(m in seq for m in motifs)


def check(seq, gc_range=(0.25, 0.75), max_run=3, motifs=()):
    """Measure a sequence and flag which constraints it satisfies."""
    gc = gc_fraction(seq)
    run = max_homopolymer_run(seq)
    return {
        "gc_fraction": gc,
        "max_homopolymer_run": run,
        "gc_ok": gc_range[0] <= gc <= gc_range[1],
        "homopolymer_ok": run <= max_run,
        "motif_ok": not has_forbidden_motif(seq, motifs),
    }


if __name__ == "__main__":
    demo = "ACGTACGTGGGCAT"
    print(demo, check(demo))
