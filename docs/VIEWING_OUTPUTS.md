# Viewing generated plots

The estimator writes its figure to `PREFIX_plot.png`. Opening that file is platform-specific; the analysis itself does not require a graphical desktop.

```bash
# macOS
open my_species_plot.png

# Linux desktop
xdg-open my_species_plot.png

# Windows PowerShell
start my_species_plot.png
```

On a headless server or HPC login node, do not expect any of these commands to display a window. Copy the PNG to a local machine (for example with `scp`) or inspect the numeric results in `PREFIX_summary.txt` and `PREFIX_mean_accumulation.csv`.
