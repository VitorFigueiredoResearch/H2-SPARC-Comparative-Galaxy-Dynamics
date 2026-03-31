# SPARC Galaxy Parameter Extraction for H1

This directory contains a self-contained, reproducible pipeline that generates
the baryonic parameter table used in the H1 Baryon-Convolved Potential Model.

## Purpose

This tool is a reproducibility utility for the **upstream baryonic input stage**
of the H1 pipeline. It converts publicly available SPARC photometric tables into
the structural parameter file (`galaxies.csv`) consumed by the H1 model runner.

It is not a rotation-curve fitting tool and does not produce model velocities.
## Before You Use This Tool

This utility is provided as a reproducibility aid for reconstructing the
SPARC-derived H1 galaxy parameter table under the documented assumptions
of this project. It is not a general photometric extraction engine and
is not validated for arbitrary external surveys.

## Please note the following before use:

- **Cite SPARC appropriately.** The upstream source tables originate from
  the SPARC project (Lelli, McGaugh & Schombert 2016), and any scientific
  use of the extracted parameters should acknowledge that source.

- **Check source-data redistribution terms.** If you redistribute this
  package publicly, verify that the bundled SPARC/CDS source tables may
  be redistributed under the applicable data-usage terms.

- **Check assumptions before downstream use.** The generated parameters
  depend on project-specific constants, fallback logic, and SPARC-specific
  formatting assumptions documented below. Users are responsible for
  confirming that those assumptions are appropriate for their own use case.

- **Do not treat this as a full H1/H2/H3 external pipeline.** This tool
  reconstructs one upstream input layer only. It does not generate model
  velocities or replace the full internal analysis workflow.
---

## Contents

- `extract_sparc_params_v5.py` — main extraction script
- `diagnose_table1.py` — diagnostic validator; checks which `.mrt` file
  contains real SPARC galaxy data and verifies that the L[3.6], Rd, and
  MHI byte-slice offsets resolve to numeric values before a full run
- `RAW/` — unmodified SPARC source data:
    * `table1.mrt` — SPARC Galaxy Sample (Lelli+2016); **required**
    * `table2.mrt` — rotation curve data (present for completeness)
    * `Bulges.mrt` — bulge luminosities; optional; used as Mstar fallback
    * `wise_ii table1.mrt` — WISE II stellar mass estimates; optional; note
      the space in the filename is present in the source archive
    * `rotmod/` — 175 HI rotation curve `.dat` files; optional; used as
      Rd fallback if Rd is absent from table1.mrt for a given galaxy
- `output/` — reference outputs committed to repository:
    * `galaxies.csv` — baryonic parameter table (name, Rd_star, Mstar,
      hz_star, Rd_gas, Mgas, hz_gas); 175-row reference file
    * `galaxies_h1_lines.JSON` — same data in H1 JSON-lines format
      (note: file extension is uppercase `.JSON`)
    * `missing_report.txt` — lists galaxies where any field defaulted to
      0.0; **generated at runtime only**, not committed to the repository

---

## How to Run

> **Important:** `extract_sparc_params_v5.py` uses bare relative paths and
> expects all `.mrt` source files and the `rotmod/` folder to exist in the
> **current working directory** when the script is invoked.
> The source files reside in `RAW/`, so the script must be run from there:

```bash
cd tools/sparc_extractor/RAW
python ../extract_sparc_params_v5.py
```

Output files (`galaxies.csv`, `galaxies_h1_lines.json`, `missing_report.txt`)
will be written into the directory from which the script is run (i.e., `RAW/`).
Move them to `output/` if needed to match the committed reference structure.

**Note on output filename case:**
The script writes `galaxies_h1_lines.json` (lowercase extension). The committed
reference file in `output/` is `galaxies_h1_lines.JSON` (uppercase). On
case-sensitive filesystems (Linux) these are distinct filenames. On
case-insensitive filesystems (Windows, macOS) they resolve to the same file.
When verifying output, compare file contents, not filenames.

**Verification run:**

```bash
cd tools/sparc_extractor/RAW
python ../diagnose_table1.py
# Review table1_diagnostic.log to confirm byte-slice parsing before extraction
python ../extract_sparc_params_v5.py
# Compare galaxies.csv against committed output/galaxies.csv
```

---

## Fallback Logic and Assumptions

The extractor applies the following H1 mass-to-light constants (hardcoded;
do not modify without understanding the downstream H1 parameter implications):

| Constant | Value | Role |
|----------|-------|------|
| `ML_DISK` | 0.5 | Disk mass-to-light ratio |
| `ML_BULGE` | 0.7 | Bulge mass-to-light ratio |
| `HELIUM` | 1.33 | He correction factor for gas mass |
| `HZ_STAR_FACTOR` | 0.2 | hz_star = 0.2 × Rd |
| `RD_GAS_FACTOR` | 1.8 | Rd_gas = 1.8 × Rd |
| `HZ_GAS` | 0.15 | Fixed gas scale height [kpc] |

**Stellar mass fallback chain:**
1. L[3.6] from `table1.mrt` (primary)
2. WISE II log(Mstar) from `wise_ii table1.mrt` → converted to implied L[3.6]
   (treated as disk-like; bulge fraction unknown)
3. Default: Mstar = 0.0 (flagged in `missing_report.txt`)

**Rd (disk scale radius) fallback chain:**
1. Rd from `table1.mrt` (primary)
2. Header-embedded Rd from the galaxy's `rotmod/` `.dat` file
3. Default: Rd_star = 0.0 (flagged in `missing_report.txt`)

**Bulge luminosity:**
Read from `Bulges.mrt` if present; added to Mstar via ML_BULGE.

---

## Limitations and Safety Notes

1. **Portability:** This tool is designed for SPARC photometric data
   (Lelli, McGaugh & Schombert 2016). The byte-slice offsets and column
   conventions are SPARC-specific. It is not intended for use with other
   photometric surveys or galaxy catalogues without documented adaptation.

2. **H1 constant dependency:** The mass-to-light and scale-height constants
   are tuned for the H1 Baryon-Convolved Potential framework. Using the
   output `galaxies.csv` with a different potential model requires
   independent validation of these constants.

3. **Output filename case sensitivity:** The script writes lowercase
   `galaxies_h1_lines.json`; the committed reference file has uppercase
   extension `.JSON`. See the run instructions above.

4. **Working directory sensitivity:** The script must be run from the
   directory that contains the `.mrt` files. Running from `sparc_extractor/`
   (without `cd RAW`) will fail with "Table1.mrt not found".

5. **Missing data:** Galaxies with absent L[3.6], Rd, or MHI will have
   those fields set to 0.0 in the output. Check `missing_report.txt`
   after any run and assess whether 0.0-field galaxies are acceptable
   for downstream use.

6. **Source data redistribution:** The `.mrt` files in `RAW/` originate
   from the SPARC project (CDS/VizieR). Users who wish to redistribute
   this repository publicly should verify the applicable data usage terms
   for SPARC source tables before including `RAW/` contents in a public
   release.

---

## Role of diagnose_table1.py

`diagnose_table1.py` is a pre-run diagnostic that validates the `.mrt`
file structure before a full extraction run. It checks:

- which `.mrt` file in the current directory contains real SPARC galaxy rows
- whether L[3.6], Rd, and MHI byte slices produce numeric values
- whether a zero-Mstar condition would be triggered

Run it before `extract_sparc_params_v5.py` if you have modified or
replaced any source `.mrt` files. Output is written to `table1_diagnostic.log`.

---

## Citations

If you use the H1 model or the SPARC-based galaxy parameters, please cite:

**SPARC (Lelli, McGaugh & Schombert 2016)**
Lelli, F., McGaugh, S. S., & Schombert, J. M. (2016).
*SPARC: Mass Models for 175 Disk Galaxies with Spitzer Photometry and
Accurate Rotation Curves.*
The Astronomical Journal, 152, 157.
doi: 10.3847/0004-6256/152/6/157

**H1 Model (Figueiredo 2025)**
Figueiredo, V. M. F. (2025). *Baryon-Convolved Effective Potential (H1).*
Zenodo. doi: 10.5281/zenodo.16967259
Repository: https://github.com/VitorFigueiredoResearch/Baryon-Convolved-Effective-Potential
