"""
mond_velocity.py
================
MOND/RAR velocity computation module.

Implements the Radial Acceleration Relation (RAR) interpolation function
from McGaugh, Lelli & Schombert (2016), as used by Li et al. (2020) to
fit all 175 SPARC galaxies.

Formula:
  g_obs = g_bar / (1 - exp(-sqrt(g_bar / g_dag)))

where:
  g_bar = baryonic centripetal acceleration = V_bar^2(r) / r
  g_dag = fundamental MOND acceleration scale (default 1.2e-10 m/s^2)
  g_obs = total predicted centripetal acceleration

Usage:
  from mond_velocity import mond_velocity, g_dag_SI, SI_to_kpc_kms

Run standalone to produce validation plot for NGC3198.
"""

import numpy as np
import os
import sys

# ── Physical constants ─────────────────────────────────────────────────────────
# Fundamental MOND acceleration scale (McGaugh+2016)
g_dag_SI = 1.2e-10          # m/s^2

# Unit conversion: m/s^2 → km^2 s^-2 kpc^-1
#   1 kpc = 3.085677581e19 m
#   g [m/s^2] * 3.085677581e19 [m/kpc] * (1e-3)^2 [(km/m)^2] = g [km^2/s^2/kpc]
_KPC_IN_M = 3.085677581e19   # m per kpc
_MS2_TO_KPC_KMS2 = _KPC_IN_M * 1e-6   # (km/s)^2/kpc per m/s^2

# g_dag in working units (km^2 s^-2 kpc^-1)
g_dag_kpc_kms = g_dag_SI * _MS2_TO_KPC_KMS2   # ≈ 3703.2 km^2/s^2/kpc


def rar_g_obs(g_bar_kpc, g_dag=g_dag_kpc_kms):
    """
    Apply the RAR interpolation function.

    Parameters
    ----------
    g_bar_kpc : array-like
        Baryonic centripetal acceleration in km^2 s^-2 kpc^-1.
        g_bar = V_bar^2(r) / r  where V_bar in km/s, r in kpc.
    g_dag : float, optional
        MOND acceleration scale in km^2 s^-2 kpc^-1.
        Default: 1.2e-10 m/s^2 converted to km^2/s^2/kpc.

    Returns
    -------
    g_obs : ndarray
        Total predicted acceleration in km^2 s^-2 kpc^-1.
    """
    g_bar_kpc = np.asarray(g_bar_kpc, dtype=float)
    # Avoid division by zero for r=0 points
    safe = g_bar_kpc > 0.0
    g_obs = np.zeros_like(g_bar_kpc)
    x = np.sqrt(g_bar_kpc[safe] / g_dag)
    g_obs[safe] = g_bar_kpc[safe] / (1.0 - np.exp(-x))
    return g_obs


def mond_velocity(r_kpc, V_gas_kms, V_disk_kms, V_bul_kms,
                  Ydisk=1.0, Ybul=1.0, g_dag=g_dag_kpc_kms):
    """
    Compute MOND/RAR total circular velocity.

    Parameters
    ----------
    r_kpc : array-like
        Galactocentric radii in kpc.
    V_gas_kms : array-like
        Gas component circular velocity in km/s (signed, squared correctly).
    V_disk_kms : array-like
        Stellar disk circular velocity in km/s.
    V_bul_kms : array-like
        Stellar bulge circular velocity in km/s.
    Ydisk : float
        Stellar M/L ratio for disk (scales V_disk^2).
    Ybul : float
        Stellar M/L ratio for bulge (scales V_bul^2).
    g_dag : float
        MOND acceleration scale in km^2 s^-2 kpc^-1.

    Returns
    -------
    V_tot : ndarray
        Total MOND circular velocity in km/s.
    V_bar : ndarray
        Baryonic circular velocity in km/s.
    """
    r   = np.asarray(r_kpc,      dtype=float)
    Vg  = np.asarray(V_gas_kms,  dtype=float)
    Vd  = np.asarray(V_disk_kms, dtype=float)
    Vb  = np.asarray(V_bul_kms,  dtype=float)

    # Baryonic velocity squared (preserve sign convention from SPARC rotmod)
    V2_bar = (np.sign(Vg) * Vg**2 +
              Ydisk * np.sign(Vd) * Vd**2 +
              Ybul  * np.sign(Vb) * Vb**2)
    V_bar = np.sqrt(np.maximum(V2_bar, 0.0))

    # Baryonic acceleration
    r_safe = np.where(r > 0, r, np.nan)
    g_bar  = V_bar**2 / r_safe     # km^2/s^2 / kpc

    # RAR predicted acceleration
    g_obs  = rar_g_obs(np.nan_to_num(g_bar, nan=0.0), g_dag=g_dag)

    # Total velocity
    V_tot  = np.sqrt(np.maximum(g_obs * r_safe, 0.0))
    V_tot  = np.where(np.isnan(V_tot), 0.0, V_tot)

    return V_tot, V_bar


# ── Validation ─────────────────────────────────────────────────────────────────

def _load_rotmod(rotmod_path):
    """Parse SPARC rotmod file → R, Vobs, errV, Vgas, Vdisk, Vbul (all km/s)."""
    data = []
    with open(rotmod_path, encoding='utf-8', errors='replace') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#') or line.startswith('!'):
                continue
            try:
                vals = [float(x) for x in line.split()]
                if len(vals) >= 6:
                    data.append(vals[:6])
            except ValueError:
                continue
    if not data:
        raise ValueError(f"No data in {rotmod_path}")
    return np.array(data).T   # R, Vobs, errV, Vgas, Vdisk, Vbul


def _parse_rar_table(rar_path):
    """
    Parse parameter_RAR.mrt into a dict keyed by galaxy name.
    Returns: {name: {'Ydisk': ..., 'Ybul': ..., 'Chi': ...}}
    """
    params = {}
    with open(rar_path, encoding='utf-8', errors='replace') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            name = line[:14].strip()
            vals = line[14:].split()
            if len(vals) < 9:
                continue
            try:
                params[name] = {
                    'Ydisk': float(vals[0]),
                    'e_Ydisk': float(vals[1]),
                    'Ybul':  float(vals[2]),
                    'e_Ybul': float(vals[3]),
                    'D':     float(vals[4]),
                    'inc':   float(vals[6]),
                    'Chi':   float(vals[8]) if vals[8] not in ('inf','nan','--') else np.inf,
                }
            except (ValueError, IndexError):
                continue
    return params


def validate_ngc3198(repo_root=None):
    """
    Validation for NGC3198: plot V_obs, V_bar, V_MOND_baseline,
    and perturbed V_MOND (g_dag ±20%).
    Saves figure to comparative_analysis/mond/figures/validation_NGC3198.png
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    if repo_root is None:
        repo_root = os.path.abspath(
            os.path.join(os.path.dirname(__file__), '..', '..'))

    rotmod = os.path.join(repo_root, 'data', 'sparc', 'NGC3198_rotmod.dat')
    rar_tab = os.path.join(repo_root, 'data', 'nfw', 'Fits', 'ByModel',
                           'Table', 'parameter_RAR.mrt')

    print("[validation] Loading NGC3198 rotmod...")
    R, Vobs, errV, Vgas, Vdisk, Vbul = _load_rotmod(rotmod)

    print("[validation] Loading RAR parameters...")
    params = _parse_rar_table(rar_tab)
    p = params.get('NGC3198')
    if p is None:
        print("  NGC3198 not found in RAR table — cannot validate.")
        return
    Ydisk, Ybul = p['Ydisk'], p['Ybul']
    print(f"  Ydisk={Ydisk:.2f}, Ybul={Ybul:.2f}, Chi={p['Chi']:.2f}")

    # Baseline MOND
    Vtot_base, Vbar = mond_velocity(R, Vgas, Vdisk, Vbul, Ydisk, Ybul)

    # Perturbed g_dag
    Vtot_p20, _ = mond_velocity(R, Vgas, Vdisk, Vbul, Ydisk, Ybul,
                                 g_dag=g_dag_kpc_kms * 1.20)
    Vtot_m20, _ = mond_velocity(R, Vgas, Vdisk, Vbul, Ydisk, Ybul,
                                 g_dag=g_dag_kpc_kms * 0.80)

    # Inner region
    R_max = R.max()
    inner = R < 0.5 * R_max
    sigma_base  = np.sqrt(np.mean((np.log10(Vtot_base[inner]) - np.log10(Vobs[inner]))**2))
    sigma_p20   = np.sqrt(np.mean((np.log10(Vtot_p20[inner])  - np.log10(Vobs[inner]))**2))
    sigma_m20   = np.sqrt(np.mean((np.log10(Vtot_m20[inner])  - np.log10(Vobs[inner]))**2))
    print(f"  Inner region: R < {0.5*R_max:.1f} kpc  ({inner.sum()} points)")
    print(f"  sigma_baseline  = {sigma_base:.4f} dex")
    print(f"  sigma (g†+20%)  = {sigma_p20:.4f} dex  (Δσ = {sigma_p20-sigma_base:+.4f})")
    print(f"  sigma (g†-20%)  = {sigma_m20:.4f} dex  (Δσ = {sigma_m20-sigma_base:+.4f})")

    # Plot
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.errorbar(R, Vobs, yerr=errV, fmt='ko', ms=4, capsize=2, label='$V_{obs}$', zorder=5)
    ax.plot(R, Vbar,      'g--',  lw=1.5, label='$V_{bar}$ (RAR Ydisk)')
    ax.plot(R, Vtot_base, 'b-',   lw=2,   label='$V_{MOND}$ baseline')
    ax.plot(R, Vtot_p20,  'b:',   lw=1.5, label='$V_{MOND}$ $g_\dagger$+20%')
    ax.plot(R, Vtot_m20,  'b--',  lw=1.5, label='$V_{MOND}$ $g_\dagger$−20%')
    ax.axvline(0.5 * R_max, color='gray', lw=1, ls='--', label='inner/outer boundary')
    ax.set_xlabel('R [kpc]')
    ax.set_ylabel('V [km/s]')
    ax.set_title(f'NGC3198  RAR validation\n'
                 f'$\\sigma_{{base}}$={sigma_base:.4f} dex  '
                 f'$\\Delta\\sigma_{{+20\\%}}$={sigma_p20-sigma_base:+.4f}  '
                 f'$\\Delta\\sigma_{{-20\\%}}$={sigma_m20-sigma_base:+.4f}')
    ax.legend(fontsize=8)
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    plt.tight_layout()
    out = os.path.join(repo_root, 'comparative_analysis', 'mond',
                       'figures', 'validation_NGC3198.png')
    plt.savefig(out, dpi=150)
    plt.close()
    print(f"  Validation plot saved: {out}")


if __name__ == '__main__':
    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
    print(f"MOND velocity module — validation run")
    print(f"g_dag = {g_dag_SI:.2e} m/s^2  =  {g_dag_kpc_kms:.1f} km^2/s^2/kpc")
    validate_ngc3198(repo_root)
