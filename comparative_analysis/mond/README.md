# MOND/RAR Comparative Perturbation Analysis

## What this analysis does

Tests whether bounded parameter perturbations of the MOND Radial Acceleration
Relation (RAR) model yield a comparable inner-region scatter-neutrality pattern
to the H2 adaptive kernel result.

The H2 result: bounded adaptation of a nonlocal kernel in inner (baryon-dominated)
regions produces near-zero change in the inner-region log10-RMS scatter metric.

The NFW result (previous stage): bounded NFW halo parameter perturbations do NOT
reproduce the same near-neutral pattern — NFW max|Δσ| ≈ 0.04–0.07 dex depending
on regime, substantially larger than H2.

The MOND result (this analysis): MOND/RAR bounded perturbations of g† (±10%, ±20%)
and Ydisk (±10%, ±20%) yield max|Δσ| ≈ 0.028–0.033 dex (median by regime),
substantially smaller than NFW but still larger than the H2 null result.

## Required input files

| File | Description |
|------|-------------|
| `data/nfw/Fits/ByModel/Table/parameter_RAR.mrt` | Li et al. 2020 RAR fit parameters |
| `data/sparc/{galaxy}_rotmod.dat` | SPARC rotation curves (all 175 galaxies) |
| `H2_PUBLICATION_RELEASE/fleet_expansion_80galaxy/fleet_summary_80galaxy.csv` | H2 fleet list |
| `comparative_analysis/nfw/nfw_perturbation_summary.csv` | NFW results (reference) |
| `comparative_analysis/metric_harmonization/h2_metric_summary.csv` | H2 harmonized results |

## MOND convention

- **Model:** Radial Acceleration Relation (McGaugh, Lelli & Schombert 2016)
- **Formula:** `g_obs = g_bar / (1 - exp(-sqrt(g_bar / g†)))`
- **g† (fixed):** 1.2 × 10⁻¹⁰ m/s² = 3703 km²/s²/kpc
- **Free parameters:** Ydisk, Ybul (from Li et al. 2020 RAR fits)
- **Source:** Li et al. (2020), ApJS 247, 31

## Perturbation scheme

| Parameter | Perturbation | Magnitudes |
|-----------|-------------|------------|
| g† | ±10%, ±20% around 1.2×10⁻¹⁰ m/s² | Primary |
| Ydisk | ±10%, ±20% around best-fit | Secondary |
| g† + Ydisk (simultaneous) | ±10%, ±20% | Combined |

## Inner-region definition

`R < 0.5 × R_max` — **identical** to H2 (test3_inner_scatter.py) and NFW analysis.

## Scatter metric

`sigma = RMS of [log10(V_model) - log10(V_obs)]` in dex — **identical** to H2 and NFW.

## How to run

From the repository root:

```bash
# Step 1: validate the MOND velocity module (produces validation_NGC3198.png)
python comparative_analysis/mond/mond_velocity.py

# Step 2: run the full 80-galaxy perturbation diagnostic
python comparative_analysis/mond/mond_perturbation_diagnostic.py

# Step 3: generate comparison figure and report
python comparative_analysis/mond/h2_mond_comparison.py

# Step 4: generate the diagnostic/failure report
python comparative_analysis/mond/diagnostic_report_template.py
```

## Summary of results

| Analysis | Median max|Δσ| (baryon-dom) | Median max|Δσ| (DM-dom) | Pattern |
|----------|---------------------------|------------------------|---------|
| H2 | ~0 dex (null result) | — | Near-neutral |
| MOND/RAR | 0.032 dex | 0.028 dex | Weakly sensitive |
| NFW | 0.043 dex | 0.071 dex | Clearly sensitive |

MOND bounded perturbations produce intermediate scatter sensitivity —
smaller than NFW, but not reproducing the H2 near-neutral pattern.
No galaxies exceed 0.05 dex under MOND perturbations (0% vs ~37–90% for NFW).

See `h2_mond_comparison_report.txt` for the full narrative.

## Output files

| File | Description |
|------|-------------|
| `data_inspection_report.txt` | RAR catalog format, conventions, coverage |
| `metric_and_region_report.txt` | Harmonized metric and inner-region definition |
| `perturbation_scheme_report.txt` | MOND perturbation rationale |
| `mond_velocity.py` | RAR velocity computation module |
| `mond_perturbation_diagnostic.py` | Main 80-galaxy perturbation pipeline |
| `mond_perturbation_summary.csv` | Per-perturbation results (888 rows, 74 galaxies) |
| `diagnostic_report.txt` | Failures, exclusions, assumptions |
| `h2_mond_comparison.py` | Comparison script (loads all three analyses) |
| `h2_mond_comparison_summary.csv` | Merged H2/MOND/NFW per-galaxy table |
| `h2_mond_comparison_report.txt` | Plain-language narrative report |
| `figures/validation_NGC3198.png` | MOND module validation plot |
| `figures/mond_scatter_sensitivity.png` | Main comparison figure |
