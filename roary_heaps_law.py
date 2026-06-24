#!/usr/bin/env python3
"""
roary_heaps_law.py

Estimate pangenome openness from a Roary-style gene presence/absence matrix
using Heaps' law:

    P(n) = k * n^gamma

where:
    P(n)  = cumulative number of genes observed after sampling n genomes
    k     = scaling constant
    gamma = Heaps' law exponent

Interpretation used here:
    gamma close to 0      -> nearly closed pangenome
    0 < gamma < 1         -> open pangenome, but openness varies by magnitude
    higher gamma          -> more genes continue to be discovered as genomes are added

Important:
    Different papers/software may report related parameters differently.
    Some use alpha where alpha = 1 - gamma. Always check the equation.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit


def heaps_law(n, k, gamma):
    """Heaps' law model: P(n) = k * n^gamma."""
    return k * np.power(n, gamma)


def load_presence_absence_matrix(input_file):
    """
    Load a simple Roary-like tab-separated matrix.

    Expected format:
        First column = gene name/id
        Remaining columns = genomes
        Values = 0/1 presence absence

    Example:
        Gene    Genome_A    Genome_B
        gene1   1           0
        gene2   1           1
    """
    input_file = Path(input_file)
    if not input_file.exists():
        raise FileNotFoundError(f"Input file does not exist: {input_file}")

    df = pd.read_csv(input_file, sep="\t")

    if df.shape[1] < 2:
        raise ValueError("Input file must contain one gene column and at least one genome column.")

    gene_column = df.columns[0]
    matrix_df = df.drop(columns=[gene_column])

    # Convert all values to numeric 0/1
    matrix_df = matrix_df.apply(pd.to_numeric, errors="raise")

    if not np.isin(matrix_df.values, [0, 1]).all():
        raise ValueError("Presence/absence matrix must contain only 0 and 1 values.")

    # Rows should be genomes, columns should be genes
    presence_absence_matrix = matrix_df.to_numpy(dtype=int).T

    genome_names = list(matrix_df.columns)
    gene_names = list(df[gene_column])

    return presence_absence_matrix, genome_names, gene_names


def run_permutations(presence_absence_matrix, iterations, seed):
    """
    Randomly permute genome order and compute cumulative pangenome size.

    Returns:
        raw_df:
            All accumulation values from all iterations.

        mean_df:
            Mean and standard deviation of pangenome size for each genome count.
    """
    rng = np.random.default_rng(seed)

    num_genomes, num_genes = presence_absence_matrix.shape

    records = []

    for iteration in range(1, iterations + 1):
        genome_order = rng.permutation(num_genomes)
        cumulative_gene_presence = np.zeros(num_genes, dtype=int)

        for step, genome_index in enumerate(genome_order, start=1):
            cumulative_gene_presence += presence_absence_matrix[genome_index]
            cumulative_gene_count = np.count_nonzero(cumulative_gene_presence)

            records.append({
                "iteration": iteration,
                "genomes_sampled": step,
                "pangenome_size": cumulative_gene_count
            })

    raw_df = pd.DataFrame(records)

    mean_df = (
        raw_df
        .groupby("genomes_sampled", as_index=False)
        .agg(
            mean_pangenome_size=("pangenome_size", "mean"),
            sd_pangenome_size=("pangenome_size", "std")
        )
    )

    mean_df["sd_pangenome_size"] = mean_df["sd_pangenome_size"].fillna(0)

    return raw_df, mean_df


def fit_heaps_law(mean_df):
    """
    Fit Heaps' law to the mean accumulation curve.
    """
    x = mean_df["genomes_sampled"].to_numpy(dtype=float)
    y = mean_df["mean_pangenome_size"].to_numpy(dtype=float)

    # Initial guess: k = first observed pangenome size, gamma = 0.5
    p0 = [max(y[0], 1), 0.5]

    # k must be positive. gamma is allowed from 0 to 2 for numerical stability.
    pars, covariance = curve_fit(
        heaps_law,
        xdata=x,
        ydata=y,
        p0=p0,
        bounds=([0, 0], [np.inf, 2]),
        maxfev=10000
    )

    k, gamma = pars

    fitted = heaps_law(x, k, gamma)

    ss_res = np.sum((y - fitted) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    r_squared = 1 - (ss_res / ss_tot) if ss_tot != 0 else np.nan

    return k, gamma, r_squared


def interpret_gamma(gamma):
    """
    Practical interpretation of gamma for this equation: P(n) = k*n^gamma.

    These ranges are heuristic, not universal biological laws.
    """
    if gamma < 0.05:
        return "Nearly closed pangenome: very few new genes are added as more genomes are sampled."
    elif gamma < 0.20:
        return "Mostly closed or weakly open pangenome: gene discovery slows strongly."
    elif gamma < 0.50:
        return "Moderately open pangenome: new genes continue to appear, but at a decreasing rate."
    elif gamma < 1.00:
        return "Highly open pangenome: substantial gene discovery continues with additional genomes."
    else:
        return "Very high gamma: check dataset heterogeneity, taxonomy, annotation quality, or model fit."


def save_plot(mean_df, k, gamma, output_plot):
    """
    Save pangenome accumulation plot with fitted Heaps' law curve.
    """
    import matplotlib.pyplot as plt

    x = mean_df["genomes_sampled"].to_numpy(dtype=float)
    y = mean_df["mean_pangenome_size"].to_numpy(dtype=float)
    yerr = mean_df["sd_pangenome_size"].to_numpy(dtype=float)

    x_fit = np.linspace(x.min(), x.max(), 200)
    y_fit = heaps_law(x_fit, k, gamma)

    plt.figure(figsize=(8, 5))
    plt.errorbar(x, y, yerr=yerr, fmt="o", capsize=3, label="Mean accumulation ± SD")
    plt.plot(x_fit, y_fit, label=f"Heaps fit: P(n) = {k:.3f} n^{gamma:.3f}")
    plt.xlabel("Number of genomes sampled")
    plt.ylabel("Cumulative pangenome size")
    plt.title("Pangenome accumulation curve")
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_plot, dpi=300)
    plt.close()


def main():
    parser = argparse.ArgumentParser(
        description="Estimate pangenome openness from a Roary-style presence/absence matrix using Heaps' law."
    )

    parser.add_argument(
        "input",
        help="Tab-separated gene presence/absence matrix. First column must be gene IDs; remaining columns must be genomes with 0/1 values."
    )

    parser.add_argument(
        "-i", "--iterations",
        type=int,
        default=100,
        help="Number of random genome-order permutations. Default: 100"
    )

    parser.add_argument(
        "-s", "--seed",
        type=int,
        default=42,
        help="Random seed for reproducibility. Default: 42"
    )

    parser.add_argument(
        "-o", "--output-prefix",
        default="heaps_law",
        help="Output prefix. Default: heaps_law"
    )

    parser.add_argument(
        "--no-plot",
        action="store_true",
        help="Do not generate PNG plot."
    )

    args = parser.parse_args()

    if args.iterations < 1:
        raise ValueError("Iterations must be at least 1.")

    print("Loading input matrix...")
    matrix, genome_names, gene_names = load_presence_absence_matrix(args.input)

    num_genomes, num_genes = matrix.shape

    print(f"Number of genomes: {num_genomes}")
    print(f"Number of genes: {num_genes}")
    print(f"Iterations: {args.iterations}")
    print(f"Seed: {args.seed}")

    print("\nRunning genome-order permutations...")
    raw_df, mean_df = run_permutations(matrix, args.iterations, args.seed)

    print("Fitting Heaps' law to mean accumulation curve...")
    k, gamma, r_squared = fit_heaps_law(mean_df)

    interpretation = interpret_gamma(gamma)

    raw_output = f"{args.output_prefix}_raw_accumulation.csv"
    mean_output = f"{args.output_prefix}_mean_accumulation.csv"
    summary_output = f"{args.output_prefix}_summary.txt"
    plot_output = f"{args.output_prefix}_plot.png"

    raw_df.to_csv(raw_output, index=False)
    mean_df.to_csv(mean_output, index=False)

    with open(summary_output, "w") as out:
        out.write("Roary Heaps' Law Summary\n")
        out.write("========================\n\n")
        out.write(f"Input file: {args.input}\n")
        out.write(f"Number of genomes: {num_genomes}\n")
        out.write(f"Number of genes: {num_genes}\n")
        out.write(f"Iterations: {args.iterations}\n")
        out.write(f"Seed: {args.seed}\n\n")
        out.write("Model:\n")
        out.write("P(n) = k * n^gamma\n\n")
        out.write(f"k = {k:.6f}\n")
        out.write(f"gamma = {gamma:.6f}\n")
        out.write(f"R_squared = {r_squared:.6f}\n\n")
        out.write("Interpretation:\n")
        out.write(interpretation + "\n")

    if not args.no_plot:
        save_plot(mean_df, k, gamma, plot_output)

    print("\nDone.")
    print(f"k = {k:.6f}")
    print(f"gamma = {gamma:.6f}")
    print(f"R_squared = {r_squared:.6f}")
    print(f"Interpretation: {interpretation}")
    print("\nFiles written:")
    print(f"  {raw_output}")
    print(f"  {mean_output}")
    print(f"  {summary_output}")
    if not args.no_plot:
        print(f"  {plot_output}")


if __name__ == "__main__":
    main()
