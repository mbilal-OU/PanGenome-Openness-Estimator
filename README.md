# PanGenome-Openness-Estimator

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python](https://img.shields.io/badge/python-3.8%2B-blue.svg)]()
[![CI](https://github.com/mbilal-OU/PanGenome-Openness-Estimator/actions/workflows/ci.yml/badge.svg)](https://github.com/mbilal-OU/PanGenome-Openness-Estimator/actions)
[![Version](https://img.shields.io/badge/version-2.0.0-green.svg)]()

> Estimate pangenome openness from a Roary gene presence–absence matrix using Heaps' law,
> permutation-based genome accumulation curves, and publication-ready figures.

---

## Table of Contents

- [What This Does](#what-this-does)
- [The Model](#the-model)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Step-by-Step Guide](#step-by-step-guide)
  - [Step 1 — Prepare your input](#step-1--prepare-your-input)
  - [Step 2 — Run the analysis](#step-2--run-the-analysis)
  - [Step 3 — Interpret your results](#step-3--interpret-your-results)
- [All Options](#all-options)
- [Input Formats](#input-formats)
- [Output Files](#output-files)
- [Real Dataset Workflow (Roary)](#real-dataset-workflow-roary)
- [Example Results](#example-results)
- [Gamma Interpretation Guide](#gamma-interpretation-guide)
- [Troubleshooting](#troubleshooting)
- [Citation](#citation)

---

## What This Does

When you add genomes one by one to a pangenome, the total number of unique genes
(the cumulative pangenome size) grows. The rate at which it grows tells you
whether your pangenome is **open** or **closed**.

```
Open pangenome:
    Each new genome adds many new genes.
    The pangenome keeps expanding.
    Species is highly diverse — many unique genes per strain.

Closed pangenome:
    Each new genome adds very few new genes.
    The pangenome curve flattens quickly.
    Species is conserved — strains share most genes.
```

This tool:
1. Takes your Roary output (or any 0/1 gene matrix)
2. Randomly permutes genome order hundreds of times
3. Computes the mean accumulation curve
4. Fits Heaps' law to estimate **gamma** — the openness exponent
5. Generates a publication-ready accumulation plot

---

## The Model

```
P(n) = k × n^gamma
```

| Parameter | Description |
|---|---|
| `P(n)` | Cumulative pangenome size after sampling `n` genomes |
| `k` | Scaling constant |
| `gamma` | Openness exponent — the key result |
| `R²` | Goodness of fit (closer to 1.0 = better fit) |

> **Important:** Some papers and tools report `alpha = 1 - gamma`. Always check
> the exact equation before comparing values across studies.

---

## Installation

### Option 1 — pip (recommended)

```bash
git clone https://github.com/mbilal-OU/PanGenome-Openness-Estimator.git
cd PanGenome-Openness-Estimator
pip install -r requirements.txt
```

### Option 2 — conda

```bash
conda install -c conda-forge numpy pandas scipy matplotlib
git clone https://github.com/mbilal-OU/PanGenome-Openness-Estimator.git
cd PanGenome-Openness-Estimator
```

### Verify installation

```bash
python3 roary_heaps_law.py --version
# PanGenome-Openness-Estimator 2.0.0
```

---

## Quick Start

```bash
# Test with the included toy datasets
python3 roary_heaps_law.py examples/toy_open_pangenome.tsv \
    --iterations 100 --seed 42 --output-prefix open_test

python3 roary_heaps_law.py examples/toy_closed_pangenome.tsv \
    --iterations 100 --seed 42 --output-prefix closed_test
```

**For your real Roary output:**

```bash
python3 roary_heaps_law.py gene_presence_absence.csv \
    --iterations 1000 \
    --seed 42 \
    --output-prefix my_species \
    --title "My Species Pangenome"
```

---

## Step-by-Step Guide

### Step 1 — Prepare your input

**If you have Roary output:**
Roary produces `gene_presence_absence.csv` automatically.
This script reads it directly — no conversion needed.

```bash
ls your_roary_output/
# gene_presence_absence.csv   ← use this file directly
```

**If you have a custom matrix:**
Prepare a tab-separated file with this format:

```
Gene        Genome_A    Genome_B    Genome_C
gene_001    1           1           0
gene_002    1           0           1
gene_003    0           1           1
```

Rules:
- First column = gene IDs (any names)
- Remaining columns = one per genome
- Values must be `0` (absent) or `1` (present)
- Tab-separated, not comma-separated

**Convert Roary CSV to TSV (optional, for other tools):**

```bash
python3 convert_roary_to_tsv.py \
    gene_presence_absence.csv \
    output_matrix.tsv
```

---

### Step 2 — Run the analysis

```bash
python3 roary_heaps_law.py <input_file> [options]
```

**Minimum required:**

```bash
python3 roary_heaps_law.py gene_presence_absence.csv
```

**Recommended for publication:**

```bash
python3 roary_heaps_law.py gene_presence_absence.csv \
    --iterations 1000 \
    --seed 42 \
    --output-prefix deinococcus_pangenome \
    --title "Deinococcus Pangenome Accumulation Curve"
```

**Skip the plot (faster, HPC-friendly):**

```bash
python3 roary_heaps_law.py gene_presence_absence.csv \
    --iterations 1000 \
    --no-plot \
    --output-prefix results
```

---

### Step 3 — Interpret your results

After running, check your `*_summary.txt` file:

```
PanGenome-Openness-Estimator v2.0.0
==================================================
Input file   : gene_presence_absence.csv
Genomes      : 36
Genes        : 6,842
Iterations   : 1000
Seed         : 42

Model
-----
P(n) = k * n^gamma

k            : 1245.832100
gamma        : 0.242000
R_squared    : 0.998700

Interpretation
--------------
Moderately open pangenome — new genes continue to appear at a decreasing rate.
```

The key number is **gamma**. See the [Gamma Interpretation Guide](#gamma-interpretation-guide) below.

---

## All Options

```
python3 roary_heaps_law.py [input] [options]

Positional arguments:
  input                 Input file (Roary CSV or 0/1 TSV)

Optional arguments:
  -i, --iterations N    Number of random permutations (default: 100)
                        Use 1000 for publication-quality results
  -s, --seed N          Random seed for reproducibility (default: 42)
  -o, --output-prefix   Prefix for output files (default: heaps_law)
  -t, --title TEXT      Custom plot title
  --no-plot             Skip PNG plot generation
  --version             Show version and exit
  -h, --help            Show help message
```

---

## Input Formats

### Format 1 — Roary gene_presence_absence.csv (auto-detected)

```
"Gene","Non-unique Gene name","Annotation","No. isolates",...,"Genome_1","Genome_2"
"dnaA","","Chromosomal replication initiator protein DnaA","36",...,"dnaA_1","dnaA_2"
"gyrA","","DNA gyrase subunit A","36",...,"gyrA_1","gyrA_2"
```

The script automatically detects Roary format and converts presence strings to 0/1 internally.

### Format 2 — Simple 0/1 TSV

```
Gene    Genome_A    Genome_B    Genome_C    Genome_D
gene1   1           1           0           1
gene2   1           0           1           1
gene3   0           1           1           0
```

---

## Output Files

| File | Description |
|---|---|
| `PREFIX_summary.txt` | Key results: k, gamma, R², interpretation |
| `PREFIX_plot.png` | Publication-ready accumulation curve (300 dpi) |
| `PREFIX_mean_accumulation.csv` | Mean ± SD pangenome size per genome count |
| `PREFIX_raw_accumulation.csv` | All individual permutation results |

**What to report in your paper:**
- Number of genomes
- Number of genes
- Number of iterations
- Seed
- k value
- gamma value
- R²
- The accumulation plot as a figure

---

## Real Dataset Workflow (Roary)

Complete workflow from Roary output to published result:

```bash
# Step 1 — Run Roary (if not done yet)
roary -e -n -i 95 -cd 99 -p 8 gff_files/*.gff -f roary_output/

# Step 2 — Run pangenome openness estimation
python3 roary_heaps_law.py \
    roary_output/gene_presence_absence.csv \
    --iterations 1000 \
    --seed 42 \
    --output-prefix my_species \
    --title "My Species — Pangenome Accumulation Curve"

# Step 3 — Check results
cat my_species_summary.txt

# Step 4 — View plot
open my_species_plot.png
```

**What to include in your Methods section:**

> "Pangenome openness was assessed using PanGenome-Openness-Estimator v2.0.0
> (https://github.com/mbilal-OU/PanGenome-Openness-Estimator). Genome order was
> randomly permuted 1,000 times (seed = 42) to generate a mean pangenome
> accumulation curve with standard deviation. Heaps' law [P(n) = k × n^gamma]
> was fitted to the mean curve using non-linear least squares. The estimated
> gamma value was used to classify pangenome openness."

---

## Example Results

### Open Pangenome (gamma = 0.421, R² = 0.988)

The curve continues rising steeply. Each new genome adds substantial new genes.
This is characteristic of environmentally diverse or highly recombinogenic species.

![Open Pangenome](examples/toy_open_plot.png)

### Closed Pangenome (gamma = 0.067, R² = 0.975)

The curve flattens quickly. After sampling most genomes, very few new genes are found.
This is characteristic of obligate pathogens or highly conserved species.

![Closed Pangenome](examples/toy_closed_plot.png)

---

## Gamma Interpretation Guide

These ranges are practical guidelines used across the pangenomics literature.
They are not universal biological laws — always consider your biological context.

| Gamma | Classification | Example organisms |
|---|---|---|
| < 0.05 | Nearly closed | Obligate intracellular pathogens |
| 0.05–0.20 | Mostly closed / weakly open | Host-adapted pathogens |
| 0.20–0.50 | Moderately open | Many environmental bacteria |
| 0.50–1.00 | Highly open | Soil bacteria, E. coli |
| ≥ 1.00 | Check dataset | Possible taxonomy or annotation issues |

**Comparison with published values:**

| Species | Gamma | Reference |
|---|---|---|
| *Streptococcus pneumoniae* | 0.11 | Tettelin et al. 2005 |
| *Escherichia coli* | 0.42 | Touchon et al. 2009 |
| *Bacillus cereus* group | 0.37 | Rasko et al. 2007 |
| *Prochlorococcus marinus* | 0.09 | Kettler et al. 2007 |

**Note on parameterisation differences:**
Tettelin et al. (2005) originally defined Heaps' law for pangenomics as
`n(g) = kN^(-alpha)` for *new* genes per genome. The equivalent relationship
for *cumulative* pangenome size is `P(n) = k × n^gamma` where `gamma ≈ 1 - alpha`
in that convention. Always check which equation a tool uses before comparing values.

---

## Troubleshooting

| Problem | Cause | Fix |
|---|---|---|
| `FileNotFoundError` | Wrong file path | Check path with `ls` |
| `ValueError: only 0 and 1` | Non-binary values in matrix | Check for missing values or non-standard encoding |
| `No genome columns found` | Wrong CSV format | Ensure using `gene_presence_absence.csv` not a summary file |
| Gamma = 0.000 | All genes are core genes | Matrix may be filtered too aggressively |
| R² < 0.90 | Poor model fit | Increase `--iterations`, check for mixed taxonomy |
| Plot not generated | matplotlib not installed | `pip install matplotlib` |
| Very slow | Too many iterations | Reduce `--iterations` or use `--no-plot` on HPC |

---

## Citation

If you use PanGenome-Openness-Estimator in your research please cite:

**This tool:**
> Bilal M. PanGenome-Openness-Estimator: Estimate pangenome openness from
> Roary-style gene presence–absence matrices using Heaps' law and permutation-based
> analysis. GitHub. 2026. https://github.com/mbilal-OU/PanGenome-Openness-Estimator

**Heaps' law (original):**
> Heaps HS. Information Retrieval: Computational and Theoretical Aspects.
> Academic Press. 1978.

**Application of Heaps' law to pangenomics:**
> Tettelin H, Masignani V, Cieslewicz MJ, et al. Genome analysis of multiple
> pathogenic isolates of Streptococcus agalactiae: implications for the microbial
> pan-genome. *Proceedings of the National Academy of Sciences.*
> 2005;102(39):13950–13955. doi:10.1073/pnas.0506758102

**Roary (pangenome construction):**
> Page AJ, Cummins CA, Hunt M, et al. Roary: rapid large-scale prokaryote pan
> genome analysis. *Bioinformatics.* 2015;31(22):3691–3693.
> doi:10.1093/bioinformatics/btv421

---

## Author

**Muhammad Bilal**  
Department of Biological Sciences  
Oakland University  
Rochester, Michigan, USA

---

## License

MIT License — free to use, modify, and distribute with attribution.
