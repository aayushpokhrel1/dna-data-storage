"""Marker resynchronization inner code (S6): confine indel damage.

The current stack treats an indel-corrupted oligo as an erasure: one insertion or
deletion desynchronizes the whole record (every downstream base shifts), the CRC
fails, and the outer code has to reconstruct the entire oligo. This layer inserts a
known marker every `period` bases; on decode it re-anchors on the markers, so an
indel's damage is confined to the single run it lands in instead of cascading through
the rest of the oligo.

This is resynchronization, not full indel correction: the run containing the indel
still ends up a few bases wrong (a deletion leaves a base whose value is unknown).
The point, and what the S6 measurement shows, is that markers turn a catastrophic
whole-oligo desync into a small, localized error the outer code can absorb far more
cheaply. (Full in-place correction is Davey-MacKay watermark + soft decoding, out of
scope here.)

The marker is a fixed 4-mer that adds no long homopolymer, and decoding only looks
for it near each expected position, so chance occurrences of the marker inside the
data do not mislead the resync.
"""
MARKER = "ACGT"


def encode(bases, period=16, marker=MARKER):
    """Insert `marker` after every `period` data bases."""
    parts = []
    for i in range(0, len(bases), period):
        parts.append(bases[i:i + period])
        parts.append(marker)
    return "".join(parts)


def _find_marker(seq, marker, center, window):
    """Position of `marker` nearest `center` within +/- window, or None."""
    for d in range(window + 1):
        for pos in (center - d, center + d):
            if 0 <= pos and seq[pos:pos + len(marker)] == marker:
                return pos
    return None


def decode(received, n_data, period=16, marker=MARKER, window=3):
    """Recover the `n_data` data bases by re-anchoring on the markers.

    Each run is read up to its marker; a run made short by a deletion is padded and
    one made long by an insertion is truncated, so the indel affects only that run
    and every later run stays aligned.
    """
    m = len(marker)
    out = []
    i = 0
    produced = 0
    while produced < n_data:
        run_len = min(period, n_data - produced)
        pos = _find_marker(received, marker, i + run_len, window)
        if pos is None:
            run = received[i:i + run_len]      # no marker found: best effort
            i += run_len
        else:
            run = received[i:pos]              # re-anchor on the marker
            i = pos + m
        if len(run) >= run_len:
            run = run[:run_len]                # insertion: drop the overflow
        else:
            run = run + "A" * (run_len - len(run))  # deletion: pad the gap
        out.append(run)
        produced += run_len
    return "".join(out)[:n_data]


if __name__ == "__main__":
    import random

    rng = random.Random(0)
    bases, prev = [], "A"
    for _ in range(128):
        b = rng.choice([x for x in "ACGT" if x != prev]); bases.append(b); prev = b
    bases = "".join(bases)
    enc = encode(bases, period=16)
    p = len(enc) // 2
    corrupted = enc[:p] + enc[p + 1:]                 # one deletion
    resynced = decode(corrupted, len(bases), period=16)
    naive = bases[:p] + bases[p + 1:]                 # deletion without markers
    frac = lambda a, b: sum(x == y for x, y in zip(a, b)) / min(len(a), len(b))
    print(f"data bases        : {len(bases)}")
    print(f"marker overhead   : {(len(enc) - len(bases)) / len(bases):.0%}")
    print(f"after 1 deletion  : marker resync {frac(resynced, bases):.1%} correct, "
          f"naive {frac(naive, bases):.1%} correct")
