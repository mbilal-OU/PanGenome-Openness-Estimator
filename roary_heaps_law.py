#!/usr/bin/env python3
"""
roary_heaps_law.py  —  PanGenome-Openness-Estimator
=====================================================
Estimate pangenome openness from a gene presence/absence matrix
using Heaps' law and permutation-based genome accumulation curves.

Model:
    P(n) = k * n^gamma

where:
    P(n)  = cumulative pangenome size after sampling n genomes
    k     = scaling constant
    gamma = Heaps' law openness exponent

Interpretation of gamma (this parameterisation):
    gamma < 0.05          nearly closed pangenome
    0.05 – 0.20           mostly closed / weakly open
    0.20 – 0.50           moderately open
    0.50 – 1.00           highly open
    >= 1.00               check taxonomy / annotation quality

Note: some tools report alpha = 1 - gamma. Always verify the equation.

Accepted input formats:
    1. Simple 0/1 TSV  (Gene | Genome_A | Genome_B ...)
    2. Roary gene_presence_absence.csv  (auto-detected)

Author : Muhammad Bilal
         Department of Biological Sciences, Oakland University
         Rochester, Michigan, USA
Version: 2.0.0
"""

__version__ = "2.0.0"

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit


# ── Model ─────────────────────────────────────────────────────────────────────

def heaps_law(n, k, gamma):
    """Heaps' law: P(n) = k * n^gamma."""
    return k * np.power(n, gamma)


# ── Input loading ──────────────────────────────────────────────────────────────

def _is_roary_csv(path: Path) -> bool:
    """Detect whether a file is a Roary gene_presence_absence.csv."""
    with open(path, "r", errors="replace") as fh:
        header = fh.readline()
    # Roary CSV always starts with these columns
    return header.startswith('"Gene"') or header.startswith('Gene,')


def load_roary_csv(input_file: Path):
    """
    Load a Roary gene_presence_absence.csv file.

    Roary CSV format:
        Gene, Non-unique Gene name, Annotation, No. isolates, No. sequences,
        Avg sequences per isolate, Genome Fragment, Order within Fragment,
        Accessory Fragment, Accessory Order with Fragment, QC,
        Min group size nuc, Max group size nuc, Avg group size nuc,
        Genome_1, Genome_2, ...

    The first 14 columns are metadata; genome columns follow.
    Presence = any non-empty string. Absence = empty / NaN.
    """
    df = pd.read_csv(input_file, low_memory=False)

    # First 14 columns are Roary metadata
    metadata_cols = 14
    gene_col      = df.columns[0]
    genome_cols   = df.columns[metadata_cols:]

    if len(genome_cols) < 2:
        raise ValueError(
            "Roary CSV has fewer than 2 genome columns after the metadata block. "
            "Check that you are using gene_presence_absence.csv (not a summary file)."
        )

    # Convert to 0/1: presence = non-empty cell
    matrix_df = df[genome_cols].notna() & (df[genome_cols] != "")
    matrix_df = matrix_df.astype(int)

    gene_names   = list(df[gene_col])
    genome_names = list(genome_cols)
    matrix       = matrix_df.to_numpy(dtype=int).T   # shape: (genomes, genes)

    return matrix, genome_names, gene_names


def load_simple_tsv(input_file: Path):
    """
    Load a simple tab-separated 0/1 presence/absence matrix.

    Expected format:
        Gene    Genome_A    Genome_B    Genome_C
        gene1   1           1           0
        gene2   1           0           1
    """
    df = pd.read_csv(input_file, sep="\t")

    if df.shape[1] < 2:
        raise ValueError(
            "Input TSV must have at least one gene column and one genome column."
        )

    gene_col  = df.columns[0]
    matrix_df = df.drop(columns=[gene_col])
    matrix_df = matrix_df.apply(pd.to_numeric, errors="raise")

    if not np.isin(matrix_df.values, [0, 1]).all():
        raise ValueError("Presence/absence matrix must contain only 0 and 1.")

    gene_names   = list(df[gene_col])
    genome_names = list(matrix_df.columns)
    matrix       = matrix_df.to_numpy(dtype=int).T   # shape: (genomes, genes)

    return matrix, genome_names, gene_names


def load_presence_absence_matrix(input_file):
    """Auto-detect format and load the presence/absence matrix."""
    path = Path(input_file)

    if not path.exists():
        raise FileNotFoundError(f"Input file not found: {path}")

    # Auto-detect Roary CSV
    if path.suffix.lower() == ".csv" or _is_roary_csv(path):
        print(f"[INFO] Detected Roary CSV format: {path.name}")
        return load_roary_csv(path)

    print(f"[INFO] Loading as simple 0/1 TSV: {path.name}")
    return load_simple_tsv(path)


# ── Permutation analysis ───────────────────────────────────────────────────────

def run_permutations(presence_absence_matrix, iterations, seed):
    """
    Randomly permute genome order and compute cumulative pangenome size.

    Returns
    -------
    raw_df  : DataFrame with all accumulation values from all iterations
    mean_df : DataFrame with mean ± SD pangenome size per genome count
    """
    rng = np.random.default_rng(seed)
    num_genomes, num_genes = presence_absence_matrix.shape

    records = []
    for iteration in range(1, iterations + 1):
        if iteration % max(1, iterations // 10) == 0:
            print(f"  Iteration {iteration}/{iterations} ...", flush=True)

        order = rng.permutation(num_genomes)
        cumulative = np.zeros(num_genes, dtype=int)

        for step, idx in enumerate(order, start=1):
            cumulative += presence_absence_matrix[idx]
            records.append({
                "iteration":      iteration,
                "genomes_sampled": step,
                "pangenome_size":  int(np.count_nonzero(cumulative)),
            })

    raw_df = pd.DataFrame(records)

    mean_df = (
        raw_df
        .groupby("genomes_sampled", as_index=False)
        .agg(
            mean_pangenome_size=("pangenome_size", "mean"),
            sd_pangenome_size  =("pangenome_size", "std"),
        )
    )
    mean_df["sd_pangenome_size"] = mean_df["sd_pangenome_size"].fillna(0)

    return raw_df, mean_df


# ── Model fitting ──────────────────────────────────────────────────────────────

def fit_heaps_law(mean_df):
    """Fit Heaps' law to the mean accumulation curve."""
    x = mean_df["genomes_sampled"].to_numpy(dtype=float)
    y = mean_df["mean_pangenome_size"].to_numpy(dtype=float)

    p0 = [max(y[0], 1.0), 0.5]

    pars, _ = curve_fit(
        heaps_law, x, y,
        p0=p0,
        bounds=([0, 0], [np.inf, 2]),
        maxfev=10000,
    )
    k, gamma = pars

    fitted  = heaps_law(x, k, gamma)
    ss_res  = np.sum((y - fitted) ** 2)
    ss_tot  = np.sum((y - np.mean(y)) ** 2)
    r2      = float(1 - ss_res / ss_tot) if ss_tot != 0 else float("nan")

    return k, gamma, r2


# ── Interpretation ─────────────────────────────────────────────────────────────

def interpret_gamma(gamma: float) -> str:
    """Return a plain-language interpretation of gamma."""
    if gamma < 0.05:
        return "Nearly closed pangenome — very few new genes are added as more genomes are sampled."
    elif gamma < 0.20:
        return "Mostly closed / weakly open pangenome — gene discovery slows strongly."
    elif gamma < 0.50:
        return "Moderately open pangenome — new genes continue to appear at a decreasing rate."
    elif gamma < 1.00:
        return "Highly open pangenome — substantial gene discovery continues with additional genomes."
    else:
        return "Very high gamma (>=1.0) — check dataset heterogeneity, taxonomy, or annotation quality."


# ── Plot ───────────────────────────────────────────────────────────────────────

def save_plot(mean_df, k, gamma, r2, output_plot, title=None):
    """Save pangenome accumulation plot with fitted Heaps' law curve."""
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker

    x    = mean_df["genomes_sampled"].to_numpy(dtype=float)
    y    = mean_df["mean_pangenome_size"].to_numpy(dtype=float)
    yerr = mean_df["sd_pangenome_size"].to_numpy(dtype=float)

    x_fit = np.linspace(x.min(), x.max(), 400)
    y_fit = heaps_law(x_fit, k, gamma)

    fig, ax = plt.subplots(figsize=(9, 6))

    ax.fill_between(x, y - yerr, y + yerr,
                    alpha=0.18, color="#2E75B6", label="Mean ± SD")
    ax.plot(x, y, "o", color="#2E75B6", markersize=4,
            label="Mean accumulation")
    ax.plot(x_fit, y_fit, "-", color="#C00000", linewidth=2,
            label=rf"Heaps fit: $P(n) = {k:.2f}\,n^{{{gamma:.3f}}}$   $R^2={r2:.4f}$")

    ax.set_xlabel("Number of genomes sampled", fontsize=12)
    ax.set_ylabel("Cumulative pangenome size (genes)", fontsize=12)
    ax.set_title(title or "Pangenome accumulation curve", fontsize=13, pad=10)
    ax.legend(fontsize=10)
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda v, _: f"{int(v):,}"))
    ax.grid(True, linestyle="--", alpha=0.4)

    # Gamma annotation box
    interp = interpret_gamma(gamma)
    ax.annotate(
        f"γ = {gamma:.4f}\n{interp}",
        xy=(0.03, 0.96), xycoords="axes fraction",
        fontsize=8, va="top",
        bbox=dict(boxstyle="round,pad=0.4", fc="lightyellow", ec="gray", alpha=0.8),
    )

    fig.tight_layout()
    fig.savefig(output_plot, dpi=300)
    plt.close(fig)
    print(f"  Plot saved: {output_plot}")


# ── CLI ────────────────────────────────────────────────────────────────────────

def build_parser():
    parser = argparse.ArgumentParser(
        prog="roary_heaps_law.py",
        description=(
            "Estimate pangenome openness using Heaps' law.\n"
            "Accepts Roary gene_presence_absence.csv or a simple 0/1 TSV matrix.\n\n"
            "Model: P(n) = k * n^gamma"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "input",
        help=(
            "Input file. Either:\n"
            "  (1) Roary gene_presence_absence.csv  — auto-detected\n"
            "  (2) Tab-separated 0/1 matrix (Gene column + genome columns)"
        ),
    )
    parser.add_argument(
        "-i", "--iterations",
        type=int, default=100,
        help="Number of random genome-order permutations (default: 100).",
    )
    parser.add_argument(
        "-s", "--seed",
        type=int, default=42,
        help="Random seed for reproducibility (default: 42).",
    )
    parser.add_argument(
        "-o", "--output-prefix",
        default="heaps_law",
        help="Prefix for all output files (default: heaps_law).",
    )
    parser.add_argument(
        "-t", "--title",
        default=None,
        help="Custom plot title (default: 'Pangenome accumulation curve').",
    )
    parser.add_argument(
        "--no-plot",
        action="store_true",
        help="Skip PNG plot generation.",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )

    return parser


def main():
    parser = build_parser()
    args   = parser.parse_args()

    if args.iterations < 1:
        parser.error("--iterations must be at least 1.")

    # ── Load ──────────────────────────────────────────────────────────────────
    print(f"\nPanGenome-Openness-Estimator v{__version__}")
    print("=" * 50)
    print("Loading input matrix ...")
    matrix, genome_names, gene_names = load_presence_absence_matrix(args.input)
    num_genomes, num_genes = matrix.shape

    print(f"  Genomes : {num_genomes}")
    print(f"  Genes   : {num_genes:,}")
    print(f"  Iterations: {args.iterations}  |  Seed: {args.seed}")

    # ── Permutations ──────────────────────────────────────────────────────────
    print("\nRunning permutations ...")
    raw_df, mean_df = run_permutations(matrix, args.iterations, args.seed)

    # ── Fit ───────────────────────────────────────────────────────────────────
    print("\nFitting Heaps' law ...")
    k, gamma, r2 = fit_heaps_law(mean_df)
    interp        = interpret_gamma(gamma)

    # ── Save outputs ──────────────────────────────────────────────────────────
    raw_out     = f"{args.output_prefix}_raw_accumulation.csv"
    mean_out    = f"{args.output_prefix}_mean_accumulation.csv"
    summary_out = f"{args.output_prefix}_summary.txt"
    plot_out    = f"{args.output_prefix}_plot.png"

    raw_df.to_csv(raw_out,  index=False)
    mean_df.to_csv(mean_out, index=False)

    with open(summary_out, "w") as fh:
        fh.write(f"PanGenome-Openness-Estimator v{__version__}\n")
        fh.write("=" * 50 + "\n\n")
        fh.write(f"Input file   : {args.input}\n")
        fh.write(f"Genomes      : {num_genomes}\n")
        fh.write(f"Genes        : {num_genes:,}\n")
        fh.write(f"Iterations   : {args.iterations}\n")
        fh.write(f"Seed         : {args.seed}\n\n")
        fh.write("Model\n-----\n")
        fh.write("P(n) = k * n^gamma\n\n")
        fh.write(f"k            : {k:.6f}\n")
        fh.write(f"gamma        : {gamma:.6f}\n")
        fh.write(f"R_squared    : {r2:.6f}\n\n")
        fh.write("Interpretation\n--------------\n")
        fh.write(interp + "\n")

    if not args.no_plot:
        print("\nGenerating plot ...")
        save_plot(mean_df, k, gamma, r2, plot_out, title=args.title)

    # ── Print summary ─────────────────────────────────────────────────────────
    print("\n" + "=" * 50)
    print("RESULTS")
    print("=" * 50)
    print(f"  k          : {k:.6f}")
    print(f"  gamma      : {gamma:.6f}")
    print(f"  R²         : {r2:.6f}")
    print(f"\n  {interp}")
    print("\nFiles written:")
    print(f"  {raw_out}")
    print(f"  {mean_out}")
    print(f"  {summary_out}")
    if not args.no_plot:
        print(f"  {plot_out}")
    print()


if __name__ == "__main__":
    main()
