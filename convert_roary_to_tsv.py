#!/usr/bin/env python3
"""
convert_roary_to_tsv.py
=======================
Convert a Roary gene_presence_absence.csv to a simple 0/1 TSV
that can be used with roary_heaps_law.py or other tools.

Usage:
    python3 convert_roary_to_tsv.py gene_presence_absence.csv output_matrix.tsv

Note:
    roary_heaps_law.py v2.0+ accepts Roary CSV directly.
    This converter is provided for compatibility with other tools.
"""

import sys
import argparse
from pathlib import Path
import pandas as pd


def convert(input_csv: str, output_tsv: str, verbose: bool = True) -> None:
    path = Path(input_csv)
    if not path.exists():
        raise FileNotFoundError(f"Input not found: {path}")

    if verbose:
        print(f"Loading: {path}")

    df = pd.read_csv(path, low_memory=False)

    # Roary metadata columns (first 14)
    metadata_cols = 14
    gene_col      = df.columns[0]
    genome_cols   = list(df.columns[metadata_cols:])

    if len(genome_cols) < 1:
        raise ValueError("No genome columns found after metadata block.")

    if verbose:
        print(f"  Genes   : {len(df):,}")
        print(f"  Genomes : {len(genome_cols)}")

    # Build 0/1 matrix: presence = non-empty cell
    binary = df[genome_cols].notna() & (df[genome_cols] != "")
    binary = binary.astype(int)

    # Assemble output: Gene column + binary columns
    out = pd.concat([df[[gene_col]].rename(columns={gene_col: "Gene"}), binary], axis=1)

    out.to_csv(output_tsv, sep="\t", index=False)

    if verbose:
        print(f"Saved TSV: {output_tsv}")
        presence_counts = binary.sum(axis=1)
        core  = (presence_counts == len(genome_cols)).sum()
        cloud = (presence_counts == 1).sum()
        print(f"  Core genes (100%)   : {core:,}")
        print(f"  Cloud genes (1 genome): {cloud:,}")


def main():
    parser = argparse.ArgumentParser(
        description="Convert Roary gene_presence_absence.csv to 0/1 TSV matrix."
    )
    parser.add_argument("input",  help="Roary gene_presence_absence.csv")
    parser.add_argument("output", help="Output 0/1 TSV file")
    parser.add_argument("-q", "--quiet", action="store_true", help="Suppress output")
    args = parser.parse_args()

    convert(args.input, args.output, verbose=not args.quiet)


if __name__ == "__main__":
    main()
