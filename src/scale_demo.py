"""End-to-end scale demonstration (gap-fill toward the paper).

Stores a real, KB-scale file, this codec's own source, in DNA and recovers it
exactly after a realistic synthesis/sequencing channel. Everything else in the repo
runs on small synthetic payloads; this shows the pipeline holds at thousands of
oligos with a real file.

Run: python src/scale_demo.py
"""
import glob
import os
import sys
import time

sys.path.insert(0, os.path.dirname(__file__))
import channel
import codec
import constraints
import ecc

# Realistic operating point: ~200 nt oligos, low RS overhead, cited-plausible channel.
DATA_BYTES = 32
PARITY_RECORDS = 8       # per 200-record block -> ~4% redundancy
BLOCK_RECORDS = 200
COVERAGE = 5
P_SUB = 0.002            # within the ~0.5% substitution range of synthesis/sequencing
P_DROP = 0.01            # 1% oligo dropout


def load_payload():
    """The repo's own source files, concatenated, as the file to store."""
    return b"".join(open(f, "rb").read() for f in sorted(glob.glob(
        os.path.join(os.path.dirname(__file__), "*.py"))))


def main():
    payload = load_payload()
    t0 = time.time()
    oligos, meta = ecc.encode(payload, data_bytes=DATA_BYTES,
                              parity_records=PARITY_RECORDS,
                              block_records=BLOCK_RECORDS, base=codec)
    enc_t = time.time() - t0

    total_nt = sum(len(o) for o in oligos)
    per_oligo_run = max(constraints.max_homopolymer_run(o) for o in oligos)
    gc = constraints.gc_fraction("".join(oligos))
    n_data = -(-len(payload) // DATA_BYTES)
    overhead = (len(oligos) - n_data) / n_data

    t0 = time.time()
    reads = channel.corrupt(oligos, p_sub=P_SUB, p_drop=P_DROP,
                            coverage=COVERAGE, seed=0)
    recovered = ecc.decode(reads, meta, base=codec)
    dec_t = time.time() - t0
    ok = recovered == payload

    print(f"file stored     : {len(payload):,} bytes (this codec's own source)")
    print(f"oligos          : {len(oligos):,} at ~{total_nt // len(oligos)} nt "
          f"({total_nt:,} nt total)")
    print(f"density         : {len(payload) * 8 / total_nt:.3f} bits/nt")
    print(f"redundancy      : {overhead:.1%} parity oligos + a 4-byte CRC each")
    print(f"constraints     : GC {gc:.3f}, max homopolymer run {per_oligo_run} (per oligo)")
    print(f"channel         : {P_SUB:.1%} substitution, {P_DROP:.1%} dropout, "
          f"coverage {COVERAGE} -> {len(reads):,} reads")
    print(f"recovered exact : {ok}")
    print(f"timing          : encode {enc_t:.1f}s, channel+decode {dec_t:.1f}s")
    assert ok, "scale demo failed to recover the file exactly"


if __name__ == "__main__":
    main()
