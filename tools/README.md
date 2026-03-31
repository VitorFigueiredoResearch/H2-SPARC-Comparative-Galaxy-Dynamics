# H2-SPARC Tools — Reproducibility Utilities

**Repository:** `H2-SPARC-Comparative-Galaxy-Dynamics`
**Paper:** *Inner-Region Scatter Response to Bounded Perturbations: A Comparative Analysis of 74 SPARC Galaxies*

---

## Overview

This directory contains two standalone reproducibility utilities. They address
different stages of the H2 pipeline and are designed to be used independently.
They share no code and no data formats.

```
SPARC raw data
     │
     ▼
[Tool A: SPARC Extractor]  →  galaxies.csv / galaxies_h1_lines.JSON
     │
     │  (H1 pipeline — external, not in this directory)
     ▼
H1 baseline model velocities
     │
     │  (H2 pipeline — Phase 3 + Phase 4, not in this directory)
     ▼
H1 + H2 model velocities at observed radii
     │
     ▼
[Tool B: H2 Metric Utility]  →  |Δσ| metric result
```

Neither tool runs the full H1 or H2 pipeline. They bookend it:
Tool A generates the upstream baryonic inputs; Tool B verifies the downstream metric.

---

## Tool A — SPARC Parameter Extractor

**Location:** `sparc_extractor/`

**Purpose:** Converts raw SPARC photometric tables (Lelli+2016) into the
baryonic parameter table used to initialise the H1 model. This is the
upstream reproducibility step: re-running this tool regenerates
`galaxies.csv` (Rd_star, Mstar, hz_star, Rd_gas, Mgas, hz_gas for all
175 SPARC galaxies) from publicly available CDS data.

**Primary script:** `sparc_extractor/extract_sparc_params_v5.py`

**Input:** Raw SPARC `.mrt` files in `sparc_extractor/RAW/`

**Output:** `sparc_extractor/output/galaxies.csv` and `galaxies_h1_lines.JSON`

**Quick start:**

> **Note:** The script uses relative paths. It must be run from inside `RAW/`,
> where the `.mrt` source files reside:

```bash
cd sparc_extractor/RAW
python ../extract_sparc_params_v5.py
```

Do not run from `sparc_extractor/` directly — the script will fail to find
`table1.mrt`. See `sparc_extractor/README.md` for full run instructions.

**Who should use this:** Anyone reproducing the H1 baryonic parameter table
from scratch, or verifying that the `galaxies.csv` in the main repository
matches the SPARC source data.

**Diagnostic validator:** `sparc_extractor/diagnose_table1.py` — checks
which `.mrt` file contains the real SPARC galaxy data and validates that
L[3.6], Rd, and MHI byte slices parse correctly.

See `sparc_extractor/README.md` for full documentation.

---

## Tool B — H2 Inner-Region Scatter Metric Utility

**Location:** `h2_diagnostic_tool.py` (this directory)

**Purpose:** Evaluates the manuscript's log₁₀-RMS inner-region scatter
sensitivity metric |Δσ| = |σ_H2 − σ_H1| given user-supplied model
velocities. This is the downstream reproducibility step: given pre-computed
H1 and H2 model velocities, verify the metric matches the manuscript
definition.

**Metric:**
```
σ = sqrt( mean( (log10(V_model) − log10(V_obs))² ) )
```
computed over R < 0.5 × R_max (inner 50% of the observed radial range).

**Manuscript anchor:** 74-galaxy combined median |Δσ| = **0.000478 dex**

**Quick start:**

```bash
python h2_diagnostic_tool.py --input example_input.csv
```

**API usage:**

```python
from h2_diagnostic_tool import evaluate
result = evaluate(R_kpc=[...], V_obs_kms=[...], V_H1_kms=[...], V_H2_kms=[...])
print(result['abs_delta_sigma_dex'])
```

**Input format:** CSV with columns `R_kpc`, `V_obs_kms`, `V_H1_kms`, `V_H2_kms`
(all in kpc / km/s; model velocities must be pre-interpolated onto the
observed radial grid).

**What this tool does NOT do:** Generate H2 or H1 model velocities. The full
H2 pipeline (baryonic decomposition, frozen H1 parameters, vendored H1 runner)
must be run separately to obtain V_H1 and V_H2 columns.

See the full documentation in this directory's extended README sections below,
or consult `TOOL_FEASIBILITY_REPORT.md` for a complete dependency and limitation
analysis.

---

## Scope and Release Status

| Tool | Stage | Status | Feasibility |
|------|-------|--------|-------------|
| SPARC Extractor (`sparc_extractor/`) | Upstream (baryonic inputs) | Released | CONDITIONAL GO — requires SPARC .mrt files from CDS |
| H2 Metric Utility (`h2_diagnostic_tool.py`) | Downstream (metric verification) | Released | CONDITIONAL GO — requires pre-computed H1/H2 velocities |
| Full H1 pipeline | Mid-stream | Not released | Requires vendored SPARC runner and per-galaxy frozen parameters |
| Full H2 pipeline (Phase 3 + 4) | Mid-stream | Not released | Requires V_baryon, frozen H1 params, Rd_star catalogue |

**Neither tool is a substitute for the full pipeline.** Both are bounded
utilities for reproducibility verification at specific, audited checkpoints.

---

## Files in This Directory

| File | Description |
|------|-------------|
| `h2_diagnostic_tool.py` | H2 metric evaluation utility (Tool B) |
| `example_input.csv` | NGC 3198 example input — 43 points with H1/H2 model velocities |
| `example_output.csv` | Expected output from running Tool B on `example_input.csv` |
| `TOOL_FEASIBILITY_REPORT.md` | Full dependency and feasibility audit for Tool B |
| `TOOL_INTEGRATION_REPORT.md` | Architecture decision and integration rationale for both tools |
| `sparc_extractor/` | SPARC parameter extraction bundle (Tool A) |

---

## Safety and Scope

### What these tools do

- **Tool A** reconstructs the H1 baryonic parameter table from raw SPARC
  photometric data using documented mass-to-light constants.
- **Tool B** evaluates the manuscript inner-region scatter metric (|Δσ| in dex)
  given externally supplied H1 and H2 model velocities.
- Both tools operate at specific, bounded checkpoints in the pipeline.
  They are reproducibility aids, not a general-purpose rotation-curve
  analysis framework.

### What these tools do not do

- They do not run the H1 or H2 model pipelines.
- They do not produce rotation-curve model velocities from first principles.
- They do not claim compatibility with photometric surveys other than SPARC
  without independent validation of units, assumptions, and input formats.
- They do not constitute a universal rotation-curve benchmark.

### Usage requirements

- Users must read `sparc_extractor/README.md` before running Tool A.
  The run path, fallback logic, and working-directory sensitivity are
  non-obvious and are documented there.
- Users must read the Input Format and Known Limitations sections below
  before supplying data to Tool B. Incorrect unit or column assignments
  will produce numerically plausible but meaningless output.
- Output values from Tool B are directly comparable to the manuscript
  anchor only if input velocities were produced by the identical H1 + H2
  pipeline described in the paper.

### Legacy diagnostic note

`diagnostics/test3_inner_scatter.py` in the parent repository is a legacy
internal diagnostic. It is **not** the authoritative metric path. A known
V_total column-priority bug causes that script to resolve `V_total_H1` in
place of `V_total_H2`, producing |Δσ| ≈ 0 for all galaxies. The authoritative
full-sample script is:

```
comparative_analysis/comparative_validation/compute_h2_full74.py
```

The standalone external reproducibility tool is `tools/h2_diagnostic_tool.py`.

---

## Metric Definition (Tool B Extended Reference)

```
σ = sqrt( mean( (log10(V_model) − log10(V_obs))² ) )
```

computed over the **inner region**: R < `inner_frac` × R_max

Default `inner_frac = 0.5` (inner 50% of the observed radial range).

```
|Δσ| = |σ_H2 − σ_H1|
```

Minimum inner-region points required: 3.

**Manuscript anchor:** The H2 74-galaxy combined median |Δσ| = **0.000478 dex**
(Tier A + B + C, authoritative V_total_H2 method).

---

## Input Format (Tool B)

A plain CSV file with the following **required** columns (case-sensitive):

| Column | Type | Unit | Description |
|--------|------|------|-------------|
| `R_kpc` | float | kpc | Galactocentric radius. Must be positive and unique. |
| `V_obs_kms` | float | km/s | Observed rotation-curve velocity. Must be positive. |
| `V_H1_kms` | float | km/s | H1 baseline model velocity interpolated onto the same radii. Must be positive. |
| `V_H2_kms` | float | km/s | H2 adaptive model velocity interpolated onto the same radii. Must be positive. |

**Optional column:**

| Column | Type | Unit | Description |
|--------|------|------|-------------|
| `V_err_kms` | float | km/s | Velocity uncertainty. Passed through to output but not used in metric computation. |

**Unit requirement:** All radii must be in kpc and all velocities in km/s.

**Interpolation responsibility:** If model velocities are defined on a different
radial grid from the observed data, the user must interpolate them onto the
observed grid before calling this tool. This matches the manuscript methodology.

---

## Installation

No installation required beyond standard Python packages:

```
pip install numpy pandas matplotlib
```

`matplotlib` is optional and only required if `--plot` is used.

Python ≥ 3.8 required.

---

## CLI Reference (Tool B)

```bash
# Basic metric evaluation
python h2_diagnostic_tool.py --input example_input.csv

# Save results to CSV
python h2_diagnostic_tool.py --input my_data.csv --output my_results.csv

# Save results + diagnostic figure
python h2_diagnostic_tool.py --input my_data.csv --output my_results.csv --plot

# Custom inner-region fraction (default is 0.5)
python h2_diagnostic_tool.py --input my_data.csv --inner_frac 0.4

# Verbose mode
python h2_diagnostic_tool.py --input my_data.csv --verbose
```

---

## Known Limitations (Tool B)

1. **Metric evaluation only.** Does not generate H2 or H1 model velocities.
2. **Interpolation responsibility lies with the user.** Evaluate on the exact
   radii provided; no internal interpolation is performed.
3. **Baryonic decomposition not validated for external data.** The H2 framework
   requires V_baryon (SPARC-derived); no validated procedure exists for other surveys.
4. **Not a universal benchmark.** Results are only directly comparable to
   manuscript values if input model velocities came from the identical pipeline.
5. **Sensitivity to inner-region coverage.** With N < 5 inner points the metric
   estimate is unstable. The tool enforces N ≥ 3 as a hard minimum.

---

## Version History

| Version | Date | Notes |
|---------|------|-------|
| 1.0.0 | 2026-03-29 | Initial release. Tool B only. |
| 1.1.0 | 2026-03-31 | Added Tool A (SPARC extractor) to tools directory; unified README. |
| 1.2.0 | 2026-03-31 | Safety hardening pass: corrected Tool A run path; added Safety and Scope section; documented legacy diagnostic warning; added sparc_extractor/README.md Limitations section. |
