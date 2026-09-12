"""Publication figures for the DNA-data-storage paper.

Three figures, each built only from the existing measurement functions
(`benchmark.sweep`, `benchmark.density_sweep`, `ecc.encode` +
`constraints`): no numbers are invented here. Output goes to
results/figures/ as PNGs at 150 dpi.

Headless: the Agg backend is selected before pyplot is imported, so this runs
on a machine with no display.
"""
import math
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(__file__))
import benchmark
import constraints
import ecc

FIGURES_DIR = os.path.join(os.path.dirname(__file__), os.pardir, "results", "figures")


def _save(fig, name):
    """Write `fig` to results/figures/<name> and return the path."""
    os.makedirs(FIGURES_DIR, exist_ok=True)
    path = os.path.join(FIGURES_DIR, name)
    plt.tight_layout()
    plt.savefig(path, dpi=150)
    plt.close(fig)
    return path


def fig_recovery():
    """Recovery vs per-base error rate, one line per sequencing coverage."""
    data = bytes(range(256))
    rates = [0.0, 0.005, 0.01, 0.02, 0.05, 0.1]

    fig, ax = plt.subplots()
    for cov in [1, 8, 20]:
        res = benchmark.sweep(data, rates, trials=10, seed=0,
                              ratio=benchmark.LIT_RATIO, coverage=cov,
                              data_bytes=8, parity_records=4)
        ax.plot(res["rates"], res["recovery"], marker="o", label=f"coverage {cov}")

    ax.set_xlabel("per-base error rate")
    ax.set_ylabel("fraction recovered")
    ax.set_title("Recovery vs error rate")
    ax.set_ylim(-0.02, 1.02)
    ax.legend()
    ax.grid(alpha=0.3)
    return _save(fig, "recovery_vs_rate.png")


def fig_density():
    """Information density vs payload size, against the rotating ceiling."""
    s = benchmark.density_sweep([4, 8, 16, 32, 64, 128, 256])

    fig, ax = plt.subplots()
    ax.plot(s["payload_bytes"], s["bits_per_nt"], marker="o", label="measured density")
    ax.axhline(math.log2(3), linestyle="--", color="gray",
               label="rotating ceiling (log2 3)")

    ax.set_xscale("log", base=2)
    ax.set_xlabel("payload bytes per oligo")
    ax.set_ylabel("density (bits/nt)")
    ax.set_title("Information density vs payload size")
    ax.legend()
    ax.grid(alpha=0.3)
    return _save(fig, "density.png")


def fig_constraints():
    """GC content and homopolymer-run distributions of an encoded pool."""
    oligos, _ = ecc.encode(bytes(range(256)), data_bytes=8, parity_records=4)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))

    gcs = [constraints.gc_fraction(o) for o in oligos]
    ax1.hist(gcs, bins=20)
    ax1.axvline(0.5, linestyle="--", color="gray")
    ax1.set_xlabel("GC fraction (per oligo)")
    ax1.set_ylabel("count")
    ax1.set_title("GC content")

    runs = [constraints.max_homopolymer_run(o) for o in oligos]
    # Discrete integer bars; show at least 1..3 so the (expected) single bar at
    # run == 1 reads clearly instead of one auto-scaled block spanning the axis.
    xs = list(range(1, max(max(runs), 3) + 1))
    counts = [runs.count(v) for v in xs]
    ax2.bar(xs, counts, width=0.6)
    ax2.set_xticks(xs)
    ax2.set_xlim(0.5, xs[-1] + 0.5)
    ax2.set_xlabel("max homopolymer run (per oligo)")
    ax2.set_ylabel("count")
    ax2.set_title("Homopolymer runs (rotating code guarantees 1)")

    return _save(fig, "constraints.png")


def fig_pareto():
    """Density vs robustness for the four codec x ECC families (S3 Pareto)."""
    res = benchmark.default_family_comparison()
    fam = res["families"]

    fig, ax = plt.subplots()
    # distinct marker per family so families that coincide at this operating point
    # (identical density and threshold) stay individually visible via the legend.
    markers = ["o", "s", "^", "D"]
    for (name, d), m in zip(fam.items(), markers):
        ax.scatter(d["bits_per_nt"], d["threshold"], s=110, marker=m,
                   edgecolor="black", linewidth=0.5, alpha=0.8, label=name)
    ax.set_xlabel("density (bits/nt)")
    ax.set_ylabel("robustness: max per-base rate at full recovery")
    ax.set_title(f"Density vs robustness by family "
                 f"(overhead {res['overhead']:.0f}x, coverage {res['coverage']})")
    ax.grid(alpha=0.3)
    ax.margins(0.3)
    ax.legend(title="base codec + ECC", loc="best")
    return _save(fig, "pareto.png")


def main():
    """Render all figures and verify each file was written."""
    paths = [fig_recovery(), fig_density(), fig_constraints(), fig_pareto()]
    for path in paths:
        assert os.path.exists(path), f"missing figure: {path}"
        assert os.path.getsize(path) > 0, f"empty figure: {path}"
        print(os.path.normpath(path))


if __name__ == "__main__":
    main()
