#!/usr/bin/env python3
"""
h2_diagnostic_tool.py
=====================
Standalone diagnostic utility for the H2 inner-region scatter metric.

Repository: H2-SPARC-Comparative-Galaxy-Dynamics
Paper: "Inner-Region Scatter Response to Bounded Perturbations:
        A Comparative Analysis of 74 SPARC Galaxies"

PURPOSE
-------
This tool is a reproducibility aid and diagnostic utility. It evaluates
the H2 inner-region scatter sensitivity metric (log10-RMS |Δσ| in dex)
given user-supplied observed and model rotation-curve data.

IMPORTANT SCOPE LIMITATIONS
-----------------------------
This tool does NOT:
  - generate H2 model velocities (V_H2). Those must be produced by the
    full H2 pipeline (Phase 3 + Phase 4) run on the appropriate data.
  - run the H1 baseline model. V_H1 must be externally supplied.
  - accept arbitrary survey formats without documented conversion.
  - claim compatibility with data outside the documented input format.
  - constitute a universal benchmark or proof of any underlying theory.

This tool ONLY evaluates the manuscript metric given pre-computed inputs.
Users must supply correct model velocities themselves.

METRIC DEFINITION (matches manuscript Section 3.2)
---------------------------------------------------
  σ = sqrt(mean((log10(V_model) − log10(V_obs))²))
       computed over the inner region: R < inner_frac × R_max

  |Δσ| = |σ_H2 − σ_H1|    (absolute sensitivity)

  Default inner_frac = 0.5 (inner 50% of the radial range)
  Minimum inner points required: 3

INPUT FORMAT (CSV)
------------------
Required columns (case-sensitive):
  R_kpc       float  Galactocentric radius [kpc]. Must be positive, unique, sorted ascending.
  V_obs_kms   float  Observed rotation-curve velocity [km/s]. Must be positive.
  V_H1_kms    float  Baseline H1 model velocity at the same radii [km/s]. Must be positive.
  V_H2_kms    float  Adaptive H2 model velocity at the same radii [km/s]. Must be positive.

Optional column:
  V_err_kms   float  Velocity uncertainty [km/s]. Loaded and passed through but not used
                     in the metric computation. Included for documentation purposes.

UNIT ASSUMPTIONS
----------------
  - All radii in kpc.
  - All velocities in km/s.
  - No silent unit conversions are performed.
  - If your data are in different units, convert externally before using this tool.

MANUSCRIPT REFERENCE ANCHOR
----------------------------
  H2 74-galaxy combined median |Δσ| = 0.000478 dex
  (Tier A + Tier B + Tier C, authoritative V_total_H2 method)
  This value is printed as a reference context but does NOT imply
  that your result will match it.

USAGE
-----
  python h2_diagnostic_tool.py --input example_input.csv
  python h2_diagnostic_tool.py --input my_data.csv --output my_results.csv
  python h2_diagnostic_tool.py --input my_data.csv --inner_frac 0.5 --plot --verbose

DEPENDENCIES
------------
  Python >= 3.8
  numpy, pandas
  matplotlib (optional, only required if --plot is used)
"""

from __future__ import annotations

import argparse
import sys
import os
import numpy as np
import pandas as pd

# ─────────────────────────────────────────────────────────────────────────────
# Operational constants — match manuscript definitions exactly
# ─────────────────────────────────────────────────────────────────────────────

INNER_FRAC_DEFAULT  = 0.5          # inner region: R < INNER_FRAC * R_max
MIN_INNER_POINTS    = 3            # minimum inner-region data points required
MANUSCRIPT_ANCHOR   = 0.000478     # H2 74-galaxy combined median |Δσ| (dex)
TOOL_VERSION        = "1.0.0"

REQUIRED_COLS = ["R_kpc", "V_obs_kms", "V_H1_kms", "V_H2_kms"]
OPTIONAL_COLS = ["V_err_kms"]


# ─────────────────────────────────────────────────────────────────────────────
# Input validation
# ─────────────────────────────────────────────────────────────────────────────

class InputError(ValueError):
    """Raised when input data fails validation."""
    pass


def validate_input(df: pd.DataFrame, source_label: str = "input") -> None:
    """
    Validate required columns, types, and basic physical sanity.
    Raises InputError with a clear message on any failure.
    """
    # 1. Required columns present
    missing = [c for c in REQUIRED_COLS if c not in df.columns]
    if missing:
        raise InputError(
            f"Missing required columns in {source_label}: {missing}\n"
            f"  Found columns: {list(df.columns)}\n"
            f"  Required: {REQUIRED_COLS}"
        )

    # 2. No empty dataframe
    if len(df) == 0:
        raise InputError(f"{source_label}: dataframe has zero rows.")

    # 3. Numeric types and finite values
    for col in REQUIRED_COLS:
        try:
            vals = pd.to_numeric(df[col], errors="raise")
        except Exception:
            raise InputError(f"Column '{col}' in {source_label} contains non-numeric values.")

        if vals.isna().any():
            n_nan = vals.isna().sum()
            raise InputError(
                f"Column '{col}' in {source_label} contains {n_nan} NaN/missing value(s). "
                f"All values must be finite and non-missing."
            )

        if not np.isfinite(vals.values).all():
            raise InputError(f"Column '{col}' in {source_label} contains Inf or -Inf values.")

    # 4. Physical positivity
    for col in ["R_kpc", "V_obs_kms", "V_H1_kms", "V_H2_kms"]:
        if (df[col] <= 0).any():
            n_bad = (df[col] <= 0).sum()
            raise InputError(
                f"Column '{col}' in {source_label}: {n_bad} value(s) <= 0. "
                f"All radii and velocities must be strictly positive (units: kpc and km/s)."
            )

    # 5. No duplicate radii
    if df["R_kpc"].duplicated().any():
        dupes = df["R_kpc"][df["R_kpc"].duplicated(keep=False)].unique()
        raise InputError(
            f"Duplicate values found in R_kpc in {source_label}: {dupes[:5]}. "
            f"Each radius must appear exactly once."
        )

    # 6. Minimum row count
    if len(df) < MIN_INNER_POINTS:
        raise InputError(
            f"{source_label} has only {len(df)} rows. "
            f"Minimum {MIN_INNER_POINTS} rows required."
        )

    # 7. Warn if V_err_kms is present but has issues (non-fatal)
    if "V_err_kms" in df.columns:
        try:
            errs = pd.to_numeric(df["V_err_kms"], errors="coerce")
            if errs.isna().any() or (errs <= 0).any():
                print(
                    f"  [WARNING] V_err_kms in {source_label} contains non-positive or "
                    f"missing values. This column is not used in the metric but may affect "
                    f"downstream analyses.",
                    file=sys.stderr
                )
        except Exception:
            pass


# ─────────────────────────────────────────────────────────────────────────────
# Metric computation — matches compute_h2_full74.py authoritative method
# ─────────────────────────────────────────────────────────────────────────────

def compute_rms_log10_scatter(V_model: np.ndarray, V_obs: np.ndarray) -> float:
    """
    Log10-RMS scatter between model and observed velocities.

    Formula (manuscript definition):
        σ = sqrt(mean((log10(V_model) − log10(V_obs))²))

    Parameters
    ----------
    V_model : 1D array of positive model velocities [km/s]
    V_obs   : 1D array of positive observed velocities [km/s]
               (same length as V_model)

    Returns
    -------
    sigma : float, scatter in dex
    """
    V_model = np.asarray(V_model, dtype=np.float64)
    V_obs   = np.asarray(V_obs,   dtype=np.float64)

    if len(V_model) != len(V_obs):
        raise ValueError(
            f"Length mismatch: len(V_model)={len(V_model)}, len(V_obs)={len(V_obs)}"
        )

    valid = (V_model > 0) & (V_obs > 0) & np.isfinite(V_model) & np.isfinite(V_obs)
    if valid.sum() < MIN_INNER_POINTS:
        return np.nan

    resid = np.log10(V_model[valid]) - np.log10(V_obs[valid])
    return float(np.sqrt(np.mean(resid ** 2)))


def run_metric(
    df:         pd.DataFrame,
    inner_frac: float = INNER_FRAC_DEFAULT,
) -> dict:
    """
    Compute the H2 inner-region scatter sensitivity metric.

    Parameters
    ----------
    df : validated DataFrame with REQUIRED_COLS present
    inner_frac : float, inner region defined as R < inner_frac * R_max

    Returns
    -------
    result : dict with all computed quantities
    """
    # Sort ascending by radius
    df = df.sort_values("R_kpc").reset_index(drop=True)

    R     = df["R_kpc"].values.astype(np.float64)
    V_obs = df["V_obs_kms"].values.astype(np.float64)
    V_H1  = df["V_H1_kms"].values.astype(np.float64)
    V_H2  = df["V_H2_kms"].values.astype(np.float64)

    R_max  = float(R.max())
    R_cut  = inner_frac * R_max
    inner  = R < R_cut
    n_inner = int(inner.sum())
    n_total = len(df)

    # Inner-region check
    if n_inner < MIN_INNER_POINTS:
        return {
            "status":              "FAIL",
            "status_reason":       (
                f"Only {n_inner} point(s) in inner region "
                f"(R < {R_cut:.3f} kpc = {inner_frac:.2f} × R_max). "
                f"Minimum required: {MIN_INNER_POINTS}. "
                f"Consider increasing inner_frac or providing higher-resolution data."
            ),
            "n_total":             n_total,
            "n_inner":             n_inner,
            "R_max_kpc":           R_max,
            "R_cut_kpc":           R_cut,
            "inner_frac":          inner_frac,
            "sigma_H1_dex":        np.nan,
            "sigma_H2_dex":        np.nan,
            "delta_sigma_dex":     np.nan,
            "abs_delta_sigma_dex": np.nan,
            "manuscript_anchor_74gal_median_dex": MANUSCRIPT_ANCHOR,
        }

    sigma_H1 = compute_rms_log10_scatter(V_H1[inner], V_obs[inner])
    sigma_H2 = compute_rms_log10_scatter(V_H2[inner], V_obs[inner])

    if np.isnan(sigma_H1) or np.isnan(sigma_H2):
        return {
            "status":              "FAIL",
            "status_reason":       (
                "NaN result in inner-region scatter computation. Check that V_H1 "
                "and V_H2 are both positive and finite in the inner region."
            ),
            "n_total":             n_total,
            "n_inner":             n_inner,
            "R_max_kpc":           R_max,
            "R_cut_kpc":           R_cut,
            "inner_frac":          inner_frac,
            "sigma_H1_dex":        float(sigma_H1),
            "sigma_H2_dex":        float(sigma_H2),
            "delta_sigma_dex":     np.nan,
            "abs_delta_sigma_dex": np.nan,
            "manuscript_anchor_74gal_median_dex": MANUSCRIPT_ANCHOR,
        }

    delta_sigma     = float(sigma_H2 - sigma_H1)
    abs_delta_sigma = float(abs(delta_sigma))

    return {
        "status":              "OK",
        "status_reason":       "Metric computed successfully.",
        "n_total":             n_total,
        "n_inner":             n_inner,
        "R_max_kpc":           float(R_max),
        "R_cut_kpc":           float(R_cut),
        "inner_frac":          float(inner_frac),
        "sigma_H1_dex":        sigma_H1,
        "sigma_H2_dex":        sigma_H2,
        "delta_sigma_dex":     delta_sigma,
        "abs_delta_sigma_dex": abs_delta_sigma,
        "manuscript_anchor_74gal_median_dex": MANUSCRIPT_ANCHOR,
    }


# ─────────────────────────────────────────────────────────────────────────────
# Output and reporting
# ─────────────────────────────────────────────────────────────────────────────

def print_result(result: dict, label: str = "") -> None:
    """Print a formatted summary of the metric result."""
    prefix = f"[{label}] " if label else ""
    print()
    print(f"{prefix}H2 Inner-Region Scatter Diagnostic — v{TOOL_VERSION}")
    print("-" * 60)
    print(f"  Status          : {result['status']}")
    if result["status"] != "OK":
        print(f"  Reason          : {result['status_reason']}")

    print(f"  N total points  : {result['n_total']}")
    print(f"  N inner points  : {result['n_inner']}")
    print(f"  R_max           : {result['R_max_kpc']:.3f} kpc")
    print(f"  R_cut           : {result['R_cut_kpc']:.3f} kpc  ({result['inner_frac']:.2f} × R_max)")

    if result["status"] == "OK":
        print()
        print(f"  σ_H1 (baseline) : {result['sigma_H1_dex']:.6f} dex")
        print(f"  σ_H2 (adaptive) : {result['sigma_H2_dex']:.6f} dex")
        print(f"  Δσ = σ_H2−σ_H1  : {result['delta_sigma_dex']:+.6f} dex")
        print(f"  |Δσ|            : {result['abs_delta_sigma_dex']:.6f} dex")
        print()
        print(f"  Manuscript anchor (74-gal median): {MANUSCRIPT_ANCHOR:.6f} dex")
        print()
        # Contextual interpretation
        ads = result["abs_delta_sigma_dex"]
        if ads < 0.001:
            ctx = "near-neutral (|Δσ| < 0.001 dex; H2 and H1 essentially equivalent on this dataset)"
        elif ads < 0.01:
            ctx = "sub-milli-dex range (consistent with the manuscript sample range)"
        else:
            ctx = "elevated (|Δσ| >= 0.01 dex; inspect model quality and inner region coverage)"
        print(f"  Context         : {ctx}")

    print("-" * 60)


def save_output(result: dict, df_input: pd.DataFrame, output_path: str) -> None:
    """Save metric result as a single-row CSV alongside input metadata."""
    row = {k: [v] for k, v in result.items()}
    out = pd.DataFrame(row)
    out.to_csv(output_path, index=False)
    print(f"  Output saved to : {output_path}")


def make_diagnostic_plot(df: pd.DataFrame, result: dict, plot_path: str) -> None:
    """
    Generate a two-panel diagnostic figure:
      Left:  Rotation curves — V_obs, V_H1, V_H2 vs R_kpc, with inner-region shading
      Right: log10 residuals in the inner region for H1 and H2
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("  [WARNING] matplotlib not available. Skipping plot.", file=sys.stderr)
        return

    df = df.sort_values("R_kpc").reset_index(drop=True)
    R     = df["R_kpc"].values
    V_obs = df["V_obs_kms"].values
    V_H1  = df["V_H1_kms"].values
    V_H2  = df["V_H2_kms"].values

    R_cut   = result["R_cut_kpc"]
    inner   = R < R_cut
    n_inner = result["n_inner"]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))
    fig.suptitle(
        f"H2 Diagnostic — inner-region scatter sensitivity\n"
        f"R_cut = {R_cut:.2f} kpc  ({result['inner_frac']:.2f}×R_max)  "
        f"N_inner = {n_inner}",
        fontsize=10
    )

    # Left panel: rotation curves
    ax1.axvspan(R.min(), R_cut, alpha=0.10, color="steelblue",
                label=f"Inner region (R < {R_cut:.1f} kpc)")
    ax1.plot(R, V_obs, "ko", ms=4, alpha=0.85, label="V_obs (observed)", zorder=4)
    ax1.plot(R, V_H1,  "b--", lw=1.6, alpha=0.80, label=f"V_H1 (baseline)  σ={result['sigma_H1_dex']:.4f} dex")
    ax1.plot(R, V_H2,  "g-",  lw=1.6, alpha=0.80, label=f"V_H2 (adaptive)  σ={result['sigma_H2_dex']:.4f} dex")
    ax1.set_xlabel("R [kpc]")
    ax1.set_ylabel("V [km/s]")
    ax1.set_title("Rotation curves", fontsize=9)
    ax1.legend(fontsize=8, loc="best")
    ax1.set_xlim(left=0)
    ax1.set_ylim(bottom=0)

    # Right panel: log10 residuals in inner region
    if inner.sum() >= MIN_INNER_POINTS and result["status"] == "OK":
        r_in = R[inner]
        res_H1 = np.log10(V_H1[inner]) - np.log10(V_obs[inner])
        res_H2 = np.log10(V_H2[inner]) - np.log10(V_obs[inner])

        ax2.axhline(0, color="gray", lw=1.0, ls="--", zorder=1)
        ax2.scatter(r_in, res_H1, c="royalblue", s=28, alpha=0.85, label="H1 residuals", zorder=3)
        ax2.scatter(r_in, res_H2, c="forestgreen", s=28, marker="s", alpha=0.85,
                    label="H2 residuals", zorder=3)
        ax2.axhline(result["sigma_H1_dex"],  color="royalblue",   lw=1.2, ls=":",
                    label=f"±σ_H1 = {result['sigma_H1_dex']:.4f} dex")
        ax2.axhline(-result["sigma_H1_dex"], color="royalblue",   lw=1.2, ls=":")
        ax2.axhline(result["sigma_H2_dex"],  color="forestgreen", lw=1.2, ls=":",
                    label=f"±σ_H2 = {result['sigma_H2_dex']:.4f} dex")
        ax2.axhline(-result["sigma_H2_dex"], color="forestgreen", lw=1.2, ls=":")
        ax2.set_xlabel("R [kpc]")
        ax2.set_ylabel("log10(V_model / V_obs)")
        ax2.set_title(f"|Δσ| = {result['abs_delta_sigma_dex']:.6f} dex", fontsize=9)
        ax2.legend(fontsize=8, loc="best")
    else:
        ax2.text(0.5, 0.5, "Insufficient inner-region points\nfor residual plot",
                 ha="center", va="center", transform=ax2.transAxes)

    plt.tight_layout()
    plt.savefig(plot_path, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Diagnostic plot : {plot_path}")


# ─────────────────────────────────────────────────────────────────────────────
# Command-line interface
# ─────────────────────────────────────────────────────────────────────────────

def main(argv=None) -> int:
    """
    Entry point. Returns 0 on success, 1 on error.
    """
    parser = argparse.ArgumentParser(
        description=(
            "H2 inner-region scatter diagnostic utility (v%(version)s).\n"
            "Evaluates |Δσ| = |σ_H2 − σ_H1| given user-supplied model velocities.\n"
            "Does NOT generate model velocities. See --help for input format."
        ) % {"version": TOOL_VERSION},
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--input", "-i", required=True,
        help="Path to input CSV. Required columns: R_kpc, V_obs_kms, V_H1_kms, V_H2_kms."
    )
    parser.add_argument(
        "--output", "-o", default=None,
        help="Optional path for output CSV containing metric result."
    )
    parser.add_argument(
        "--inner_frac", type=float, default=INNER_FRAC_DEFAULT,
        help=(
            f"Fraction of R_max defining the inner region. "
            f"Default: {INNER_FRAC_DEFAULT} (matches manuscript definition). "
            f"Must be in (0, 1)."
        )
    )
    parser.add_argument(
        "--plot", action="store_true",
        help="Save a diagnostic figure alongside the output. Requires matplotlib."
    )
    parser.add_argument(
        "--verbose", action="store_true",
        help="Print additional diagnostic information."
    )

    args = parser.parse_args(argv)

    # ── Validate inner_frac ────────────────────────────────────────────────
    if not (0.0 < args.inner_frac < 1.0):
        print(
            f"ERROR: --inner_frac must be in (0, 1). Got: {args.inner_frac}",
            file=sys.stderr
        )
        return 1

    # ── Load input ─────────────────────────────────────────────────────────
    if not os.path.isfile(args.input):
        print(f"ERROR: Input file not found: {args.input}", file=sys.stderr)
        return 1

    try:
        df = pd.read_csv(args.input)
    except Exception as e:
        print(f"ERROR: Could not read input CSV '{args.input}': {e}", file=sys.stderr)
        return 1

    if args.verbose:
        print(f"  Loaded {len(df)} rows from '{args.input}'")
        print(f"  Columns found: {list(df.columns)}")

    # ── Validate ───────────────────────────────────────────────────────────
    try:
        validate_input(df, source_label=os.path.basename(args.input))
    except InputError as e:
        print(f"\nVALIDATION ERROR:\n  {e}", file=sys.stderr)
        return 1

    # ── Run metric ─────────────────────────────────────────────────────────
    result = run_metric(df, inner_frac=args.inner_frac)

    # ── Print result ───────────────────────────────────────────────────────
    print_result(result, label=os.path.basename(args.input))

    # ── Save output ────────────────────────────────────────────────────────
    if args.output:
        try:
            save_output(result, df, args.output)
        except Exception as e:
            print(f"ERROR: Could not write output file '{args.output}': {e}", file=sys.stderr)
            return 1

    # ── Plot ───────────────────────────────────────────────────────────────
    if args.plot:
        plot_path = args.output.replace(".csv", "_diagnostic.png") if args.output else \
                    args.input.replace(".csv", "_diagnostic.png")
        try:
            make_diagnostic_plot(df, result, plot_path)
        except Exception as e:
            print(f"  [WARNING] Plot generation failed: {e}", file=sys.stderr)

    return 0 if result["status"] == "OK" else 1


# ─────────────────────────────────────────────────────────────────────────────
# Importable API (for use as a library module)
# ─────────────────────────────────────────────────────────────────────────────

def evaluate(
    R_kpc:      "array-like",
    V_obs_kms:  "array-like",
    V_H1_kms:   "array-like",
    V_H2_kms:   "array-like",
    inner_frac: float = INNER_FRAC_DEFAULT,
) -> dict:
    """
    Importable API: evaluate H2 scatter metric from arrays.

    Parameters
    ----------
    R_kpc, V_obs_kms, V_H1_kms, V_H2_kms : array-like
        Equal-length 1D sequences. Units: kpc and km/s. All must be positive.
    inner_frac : float
        Fraction of R_max defining inner region. Default 0.5.

    Returns
    -------
    dict with keys: status, sigma_H1_dex, sigma_H2_dex, delta_sigma_dex,
                    abs_delta_sigma_dex, n_inner, R_cut_kpc, ...
    """
    df = pd.DataFrame({
        "R_kpc":      np.asarray(R_kpc,     dtype=float),
        "V_obs_kms":  np.asarray(V_obs_kms, dtype=float),
        "V_H1_kms":   np.asarray(V_H1_kms,  dtype=float),
        "V_H2_kms":   np.asarray(V_H2_kms,  dtype=float),
    })
    validate_input(df, source_label="evaluate() call")
    return run_metric(df, inner_frac=inner_frac)


if __name__ == "__main__":
    sys.exit(main())
