# Roary Heaps' Law Pangenome Openness Estimator

This project estimates whether a pangenome is open or closed using a Roary-style gene presence/absence matrix and Heaps' law.

## Biological idea

When genomes are added one by one, the cumulative number of observed genes increases.

If new genomes keep adding many new genes, the pangenome is considered open.

If new genomes add very few new genes, the pangenome is considered closed or nearly closed.

## Model

This script fits:

```text
P(n) = k × n^gamma
```

where:

- `P(n)` = cumulative pangenome size after sampling `n` genomes
- `k` = scaling constant
- `gamma` = pangenome openness exponent

## Practical interpretation of gamma

These ranges are practical guidelines, not universal biological laws.

| Gamma value | Interpretation |
|---:|---|
| `gamma < 0.05` | Nearly closed pangenome |
| `0.05–0.20` | Mostly closed / weakly open |
| `0.20–0.50` | Moderately open |
| `0.50–1.00` | Highly open |
| `>= 1.00` | Check taxonomy, annotation, or model fit |

Important: some papers use a different parameterization, such as `alpha = 1 - gamma`. Always check the exact equation before comparing values.

## Input format

The input should be a tab-separated gene presence/absence matrix:

```text
Gene    Genome_A    Genome_B    Genome_C
gene1   1           1           0
gene2   1           0           1
gene3   0           1           1
```

The first column contains gene IDs. Remaining columns are genomes. Values must be `0` or `1`.

## Installation

```bash
git clone https://github.com/YOUR_USERNAME/roary-heaps-law.git
cd roary-heaps-law
pip install -r requirements.txt
```

## Basic usage

```bash
python roary_heaps_law.py toy_closed_pangenome.tsv --iterations 100 --seed 42 --output-prefix closed_test
```

## Run open pangenome toy example

```bash
python roary_heaps_law.py toy_open_pangenome.tsv --iterations 100 --seed 42 --output-prefix open_test
```

## Outputs

For each run, the script produces:

```text
PREFIX_raw_accumulation.csv
PREFIX_mean_accumulation.csv
PREFIX_summary.txt
PREFIX_plot.png
```

## Example closed-like toy result

In the closed-like toy dataset, most genes are shared among genomes. New genomes add only a few new genes.

Expected behavior:

```text
low gamma
curve starts to flatten
closed or weakly open interpretation
```

## Example open-like toy result

In the open-like toy dataset, each genome has several genome-specific genes. New genomes keep adding genes.

Expected behavior:

```text
higher gamma
curve continues rising
open pangenome interpretation
```

## Recommended real dataset workflow

For a real Roary output:

1. Generate gene presence/absence matrix from Roary.
2. Convert presence/absence to `0/1` if needed.
3. Run:

```bash
python roary_heaps_law.py your_matrix.tsv --iterations 1000 --seed 42 --output-prefix your_species
```

4. Report:
   - number of genomes
   - number of genes
   - number of iterations
   - `k`
   - `gamma`
   - R²
   - accumulation plot

## Citation note

This is a simple educational/research utility. For publication, compare results with established tools such as PanGP, BPGA, GET_HOMOLOGUES, or micropan.
