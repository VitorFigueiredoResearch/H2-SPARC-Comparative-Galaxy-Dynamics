# scripts/figures — H2 Figure Generation Scripts

This directory contains standalone Python scripts that generate publication
figures for the H2 adaptive nonlocal kernel paper.

---

## `generate_acceleration_space_plot.py`

### Purpose

Generates `acceleration_space_80galaxy.pdf`: an acceleration-space diagnostic
figure showing all radial data points from the 80-galaxy H2 validation sample
in `g_tot` vs `g_bar` space, overlaid on the McGaugh et al. (2016) radial
acceleration relation (RAR).

### What the figure shows

| Element | Description |
|---|---|
| Data points | Each point = one radial measurement from a SPARC rotmod file |
| Red points | Baryon-dominated galaxies (Vobs/Vb < 0.8, 39 galaxies) |
| Blue points | Balanced galaxies (0.8 ≤ Vobs/Vb < 1.5, 30 galaxies) |
| Green points | DM-dominated galaxies (Vobs/Vb ≥ 1.5, 11 galaxies) |
| Black curve | McGaugh et al. (2016) RAR: g†=1.2×10⁻¹⁰ m/s² |
| Dotted line | 1:1 Newtonian limit (g_tot = g_bar) |
| Dashed line | Characteristic acceleration g† |
| Annotation | H2 null result: \|Δσ\| < 10⁻¹² km/s for 74 measurable systems |

### Acceleration definitions

```
g_bar(r) = [V_gas(r)² + V_disk(r)² + V_bul(r)²] / r   [m/s²]
g_tot(r) = V_obs(r)² / r                                [m/s²]
```

Velocities from SPARC `*_rotmod.dat` files (columns: Rad, Vobs, errV, Vgas,
Vdisk, Vbul). Quality cut: Vobs/errV > 2 and all quantities positive/finite.

### McGaugh+16 RAR formula

```
g_tot = g_bar / (1 - exp(-sqrt(g_bar / g†)))
g† = 1.2 × 10⁻¹⁰ m/s²
```

### Data inputs

| File | Location |
|---|---|
| `fleet_summary_80galaxy.csv` | `H2_PUBLICATION_RELEASE/fleet_expansion_80galaxy/` |
| `*_rotmod.dat` (80 files) | `data/sparc/` |

### Usage

Run from any working directory (paths auto-resolved from script location):

```bash
# Default output: scripts/figures/acceleration_space_80galaxy.pdf
PYTHONUTF8=1 python scripts/figures/generate_acceleration_space_plot.py

# Custom output path
PYTHONUTF8=1 python scripts/figures/generate_acceleration_space_plot.py \
    --output path/to/output.pdf
```

### Output

`acceleration_space_80galaxy.pdf` — 6.5 × 5.5 inch, 300 dpi, vector PDF
suitable for journal submission.

### Dependencies

- Python ≥ 3.9
- numpy, matplotlib (standard scientific stack; see `requirements.txt`)

### Citation

If using this figure or script, cite:

> Figueiredo, V. M. F. (2026). *Structural Limits of State-Driven Adaptive
> Nonlocal Kernel Modulation in Galaxy Rotation Curves: A Null-Response
> Diagnostic from 80 SPARC Galaxies.* MNRAS (submitted).

> McGaugh, S. S., Lelli, F., & Schombert, J. M. (2016). *Radial Acceleration
> Relation in Rotationally Supported Galaxies.* PRL, 117, 201101.

> Lelli, F., McGaugh, S. S., & Schombert, J. M. (2016). *SPARC: Mass Models
> for 175 Disk Galaxies.* AJ, 152, 157.
