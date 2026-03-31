# scripts/Figures — Figure Generation Scripts

This directory contains the figure-generation script for the main
comparative figure used in the manuscript:

*"Inner-Region Scatter Response to Bounded Perturbations:  
A Comparative Analysis of 74 SPARC Galaxies"*
(Figueiredo, V. M. F., 2026, MNRAS submitted)

---

## `generate_three_way_comparison.py` (v3.0)

### Purpose

Generates `three_way_comparison.png`: the main comparative figure for the
MNRAS manuscript, showing inner-region scatter sensitivity across the
three tested implementations on the common 74-galaxy comparison sample.

### What the figure shows

| Panel | Content |
|---|---|
| (a) Scatter plot | Maximum absolute scatter change $\max|\Delta\sigma|$ (dex) as a function of the inner-region baryonic-dominance proxy |
| (b) Box plots | Distribution summaries for NFW, adopted RAR, and H2 on the common $N=74$ basis |

### Key features

- NFW, adopted RAR, and H2 shown on a symmetric **74-galaxy** comparison basis
- NFW median $\max|\Delta\sigma| = 0.0548$ dex
- adopted RAR median $\max|\Delta\sigma| = 0.0321$ dex
- H2 median $|\Delta\sigma| = 0.000478$ dex
- NFW Spearman correlation annotated:
  $\rho = -0.899$, $p = 1.8\times10^{-27}$, $N = 74$
- adopted RAR null result annotated:
  $\rho = +0.0485$, $p = 0.681$, $N = 74$
- H2 shown with full-sample coverage on the common comparison set

### Data inputs

| File | Location |
|---|---|
| `h2_nfw_comparison_summary.csv` | `comparative_analysis/nfw/` |
| `h2_mond_comparison_summary.csv` | `comparative_analysis/mond/` |
| `h2_74galaxy_combined_metrics.csv` | `comparative_analysis/comparative_validation/` |

### Usage

From the repository root:

```bash
python scripts/Figures/generate_three_way_comparison.py
# Default output: paper/three_way_comparison.png

Optional custom output path:

python scripts/Figures/generate_three_way_comparison.py --output path/to/output.png
Output

three_way_comparison.png — main comparative figure for the manuscript,
formatted for MNRAS single-column usage.

Dependencies
Python ≥ 3.9
numpy
matplotlib
pandas
Notes
This script corresponds to the final manuscript state and supersedes
earlier internal figure-generation routes based on partial H2 archive
coverage.
The plotted H2 values reflect the full 74-galaxy common comparison
sample used in the final comparative analysis.
For manuscript interpretation, see the final figure caption and
discussion sections in the paper source.
Citation

If using this figure or script, cite:

Figueiredo, V. M. F. (2026). Inner-Region Scatter Response to Bounded
Perturbations: A Comparative Analysis of 74 SPARC Galaxies. MNRAS (submitted).

Li, P., Lelli, F., McGaugh, S. S., & Schombert, J. M. (2020). A Comprehensive
Catalog of Dark Matter Halo Models for SPARC Galaxies. ApJS, 247, 31.

Lelli, F., McGaugh, S. S., & Schombert, J. M. (2016). SPARC: Mass Models
for 175 Disk Galaxies. AJ, 152, 157.
