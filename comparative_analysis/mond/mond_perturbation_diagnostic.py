"""
mond_perturbation_diagnostic.py
================================
Bounded MOND/RAR perturbation analysis for the 80-galaxy H2 fleet.

For each galaxy with a valid RAR fit (Li et al. 2020), this script:
  1. Loads the SPARC rotation curve (baryonic components + observed)
  2. Computes the baseline RAR/MOND model velocity
  3. Applies bounded perturbations to g_dag (±10%, ±20%) and Ydisk (±10%, ±20%)
  4. Computes inner-region scatter in the harmonized log10-RMS metric for each
  5. Records delta_sigma = sigma_perturbed - sigma_baseline

Inner-region definition: R < 0.5 * R_max  (same as H2 and NFW)
Scatter metric: RMS of log10(V_model) - log10(V_obs) [dex]

Output: comparative_analysis/mond/mond_perturbation_summary.csv

Run from repository root:
  python comparative_analysis/mond/mond_perturbation_diagnostic.py
"""

import os
import sys
import numpy as np
import pandas as pd

# ── Repository root ────────────────────────────────────────────────────────────
REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
sys.path.insert(0, os.path.join(REPO_ROOT, 'comparative_analysis', 'mond'))

from mond_velocity import mond_velocity, g_dag_kpc_kms, _load_rotmod, _parse_rar_table

# ── Paths ──────────────────────────────────────────────────────────────────────
SPARC_DIR   = os.path.join(REPO_ROOT, 'data', 'sparc')
RAR_TABLE   = os.path.join(REPO_ROOT, 'data', 'nfw', 'Fits', 'ByModel',
                            'Table', 'parameter_RAR.mrt')
FLEET_CSV   = os.path.join(REPO_ROOT, 'H2_PUBLICATION_RELEASE',
                            'fleet_expansion_80galaxy', 'fleet_summary_80galaxy.csv')
OUT_CSV     = os.path.join(REPO_ROOT, 'comparative_analysis', 'mond',
                            'mond_perturbation_summary.csv')

INNER_FRAC  = 0.5        # inner region: R < INNER_FRAC * R_max
MIN_INNER   = 3          # minimum inner points for valid scatter estimate

# Perturbation magnitudes
PERTURB_MAG = [0.10, 0.20]   # fractional perturbations


def inner_scatter(V_model, V_obs, mask):
    """Log10 RMS scatter over inner region. Returns nan if < MIN_INNER points."""
    m = mask & (V_model > 0) & (V_obs > 0)
    if m.sum() < MIN_INNER:
        return np.nan
    return np.sqrt(np.mean((np.log10(V_model[m]) - np.log10(V_obs[m]))**2))


def analyze_galaxy(galaxy, regime, p, sparc_dir):
    """
    Run full perturbation analysis for one galaxy.
    Returns a list of row dicts (one per perturbation label).
    """
    # ── Load rotmod ────────────────────────────────────────────────────────────
    rotmod_path = os.path.join(sparc_dir, f'{galaxy}_rotmod.dat')
    if not os.path.exists(rotmod_path):
        # Try alternate name (some SPARC files use underscores differently)
        candidates = [f for f in os.listdir(sparc_dir)
                      if f.lower().startswith(galaxy.lower().replace('-',''))
                      and f.endswith('_rotmod.dat')]
        if candidates:
            rotmod_path = os.path.join(sparc_dir, candidates[0])
        else:
            return None, f"rotmod not found for {galaxy}"

    try:
        R, Vobs, errV, Vgas, Vdisk, Vbul = _load_rotmod(rotmod_path)
    except Exception as e:
        return None, f"rotmod load error: {e}"

    if len(R) < MIN_INNER:
        return None, f"too few rotmod points ({len(R)})"

    Ydisk_0 = p['Ydisk']
    Ybul_0  = p['Ybul']

    # ── Inner region mask ──────────────────────────────────────────────────────
    R_max  = R.max()
    inner  = R < INNER_FRAC * R_max
    n_inner = inner.sum()

    if n_inner < MIN_INNER:
        return None, f"only {n_inner} inner points (< {MIN_INNER})"

    # ── Baseline scatter ───────────────────────────────────────────────────────
    Vtot_base, Vbar_base = mond_velocity(R, Vgas, Vdisk, Vbul,
                                          Ydisk_0, Ybul_0, g_dag_kpc_kms)
    sigma_base = inner_scatter(Vtot_base, Vobs, inner)
    if np.isnan(sigma_base):
        return None, "baseline sigma is nan (V_model or V_obs zero in inner region)"

    # Median V_bar / V_MOND in inner region
    Vtot_inner = Vtot_base[inner]
    Vbar_inner = Vbar_base[inner]
    ratio_inner = np.median(np.where(Vtot_inner > 0,
                                      Vbar_inner / np.maximum(Vtot_inner, 1e-6),
                                      np.nan))

    rows = []

    def make_row(label, delta_s, variant):
        return {
            'Galaxy':             galaxy,
            'Regime':             regime,
            'mond_variant':       variant,
            'sigma_baseline':     round(sigma_base, 6),
            'n_inner_points':     int(n_inner),
            'inner_region_defined': 'Y',
            'perturbation_label': label,
            'delta_sigma':        round(delta_s, 6),
            'V_bar_over_V_MOND_inner_median': round(float(ratio_inner) if not np.isnan(ratio_inner) else 0.0, 4),
            'notes':              '',
        }

    # ── g_dag perturbations ────────────────────────────────────────────────────
    for mag in PERTURB_MAG:
        for sign, tag in [(+1, 'p'), (-1, 'm')]:
            g_new = g_dag_kpc_kms * (1.0 + sign * mag)
            V_pert, _ = mond_velocity(R, Vgas, Vdisk, Vbul, Ydisk_0, Ybul_0, g_new)
            s = inner_scatter(V_pert, Vobs, inner)
            if np.isnan(s):
                continue
            pct = int(mag * 100)
            label = f'gdag_{tag}{pct}'
            rows.append(make_row(label, s - sigma_base, 'g_dag'))

    # ── Ydisk perturbations ────────────────────────────────────────────────────
    for mag in PERTURB_MAG:
        for sign, tag in [(+1, 'p'), (-1, 'm')]:
            Yd_new = Ydisk_0 * (1.0 + sign * mag)
            Yd_new = max(Yd_new, 1e-3)   # stay physical
            Yb_new = Ybul_0  # keep bulge fixed unless scaling both
            V_pert, _ = mond_velocity(R, Vgas, Vdisk, Vbul,
                                       Yd_new, Yb_new, g_dag_kpc_kms)
            s = inner_scatter(V_pert, Vobs, inner)
            if np.isnan(s):
                continue
            pct = int(mag * 100)
            label = f'Ydisk_{tag}{pct}'
            rows.append(make_row(label, s - sigma_base, 'Ydisk'))

    # ── Combined perturbations (g_dag + Ydisk simultaneously) ─────────────────
    for mag in PERTURB_MAG:
        for sign, tag in [(+1, 'p'), (-1, 'm')]:
            g_new  = g_dag_kpc_kms * (1.0 + sign * mag)
            Yd_new = max(Ydisk_0 * (1.0 + sign * mag), 1e-3)
            V_pert, _ = mond_velocity(R, Vgas, Vdisk, Vbul,
                                       Yd_new, Ybul_0, g_new)
            s = inner_scatter(V_pert, Vobs, inner)
            if np.isnan(s):
                continue
            pct = int(mag * 100)
            label = f'both_{tag}{pct}'
            rows.append(make_row(label, s - sigma_base, 'combined'))

    return rows, None


def main():
    print("=" * 70)
    print("MOND/RAR Bounded Perturbation Diagnostic")
    print("=" * 70)

    # Load fleet
    fleet = pd.read_csv(FLEET_CSV)
    print(f"H2 fleet: {len(fleet)} galaxies")

    # Load RAR parameters
    print(f"Loading RAR parameters from {os.path.basename(RAR_TABLE)} ...")
    rar_params = _parse_rar_table(RAR_TABLE)
    print(f"  RAR entries: {len(rar_params)}")

    all_rows  = []
    succeeded = []
    failed    = []
    no_match  = []

    for _, row in fleet.iterrows():
        galaxy = row['galaxy']
        regime = row['regime']

        # Match to RAR catalog
        p = rar_params.get(galaxy)
        if p is None:
            no_match.append(galaxy)
            continue

        rows, err = analyze_galaxy(galaxy, regime, p, SPARC_DIR)
        if rows is None:
            print(f"  SKIP {galaxy}: {err}")
            failed.append((galaxy, err))
        else:
            all_rows.extend(rows)
            succeeded.append(galaxy)
            # Quick status
            perturbs   = [r['delta_sigma'] for r in rows]
            max_abs_ds = max(abs(d) for d in perturbs) if perturbs else float('nan')
            print(f"  OK   {galaxy:15s}  regime={regime:12s}  "
                  f"sigma_base={rows[0]['sigma_baseline']:.4f}  "
                  f"max|Δσ|={max_abs_ds:.4f}")

    # ── Summary CSV ────────────────────────────────────────────────────────────
    if not all_rows:
        print("\nNo valid results — check data paths.")
        return

    df = pd.DataFrame(all_rows)

    # Add max_abs_delta_sigma per galaxy (wide summary column)
    max_ds = (df.groupby('Galaxy')['delta_sigma']
                .apply(lambda x: x.abs().max())
                .rename('max_abs_delta_sigma'))
    df = df.merge(max_ds, on='Galaxy', how='left')

    df.to_csv(OUT_CSV, index=False)
    print(f"\nSaved: {OUT_CSV}")
    print(f"  Rows: {len(df)}")
    print(f"  Galaxies succeeded: {len(succeeded)}")
    print(f"  Galaxies failed:    {len(failed)}")
    print(f"  Galaxies no match:  {len(no_match)}")

    # Per-regime summary
    print("\n── Per-regime max|Δσ| summary ──")
    summary = (df.groupby(['Galaxy', 'Regime'])['delta_sigma']
                 .apply(lambda x: x.abs().max())
                 .reset_index()
                 .rename(columns={'delta_sigma': 'max_abs_ds'}))
    for regime, grp in summary.groupby('Regime'):
        med = grp['max_abs_ds'].median()
        mn  = grp['max_abs_ds'].mean()
        mx  = grp['max_abs_ds'].max()
        n   = len(grp)
        print(f"  {regime:12s}  n={n:3d}  median={med:.4f}  mean={mn:.4f}  max={mx:.4f}")

    if no_match:
        print(f"\nNo-match galaxies: {no_match}")
    if failed:
        print(f"Failed galaxies:")
        for g, e in failed:
            print(f"  {g}: {e}")


if __name__ == '__main__':
    main()
