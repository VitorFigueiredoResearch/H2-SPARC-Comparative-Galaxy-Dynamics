"""
H2 vs NFW Comparison
======================
Compares H2 inner-region scatter sensitivity against NFW perturbation
sensitivity for matched galaxies in a harmonized log-residual RMS metric
(dex), assuming the H2 operational scatter metric is the same log10-based
inner-region RMS used in its diagnostic implementation.

Inputs:
  - H2_PUBLICATION_RELEASE/fleet_expansion_80galaxy/fleet_summary_80galaxy.csv
  - comparative_analysis/nfw/nfw_perturbation_summary.csv

Outputs:
  - comparative_analysis/nfw/h2_nfw_comparison_summary.csv
  - comparative_analysis/nfw/h2_nfw_comparison_report.txt
"""
import os
import sys
import numpy as np
import pandas as pd
from scipy import stats

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(os.path.dirname(THIS_DIR))

H2_FLEET_PATH   = os.path.join(REPO_ROOT, 'H2_PUBLICATION_RELEASE',
                                 'fleet_expansion_80galaxy', 'fleet_summary_80galaxy.csv')
NFW_SUMMARY_PATH = os.path.join(THIS_DIR, 'nfw_perturbation_summary.csv')
OUT_CSV          = os.path.join(THIS_DIR, 'h2_nfw_comparison_summary.csv')
OUT_REPORT       = os.path.join(THIS_DIR, 'h2_nfw_comparison_report.txt')


def norm_name(s):
    return s.strip().upper().replace('_', '').replace(' ', '')


def classify_suppression(max_abs_ds):
    """
    Classify NFW perturbation sensitivity:
      < 0.01 dex -> strong suppression (very insensitive)
      0.01-0.05 -> moderate
      > 0.05    -> weak suppression (sensitive to perturbations)
    """
    if np.isnan(max_abs_ds):
        return 'undefined'
    if max_abs_ds < 0.01:
        return 'strong_suppression'
    elif max_abs_ds < 0.05:
        return 'moderate'
    else:
        return 'weak_suppression'


def main():
    print("=" * 60)
    print("H2 vs NFW Comparison")
    print("=" * 60)

    # Load H2 fleet
    print(f"\nLoading H2 fleet: {H2_FLEET_PATH}")
    h2_df = pd.read_csv(H2_FLEET_PATH)
    h2_df.columns = [c.strip() for c in h2_df.columns]
    h2_df['_norm'] = h2_df['galaxy'].apply(norm_name)
    print(f"  H2 fleet galaxies: {len(h2_df)}")

    # Load NFW perturbation summary
    print(f"\nLoading NFW perturbation summary: {NFW_SUMMARY_PATH}")
    if not os.path.exists(NFW_SUMMARY_PATH):
        print(f"ERROR: {NFW_SUMMARY_PATH} not found.")
        print("Please run nfw_perturbation_diagnostic.py first.")
        return

    nfw_df = pd.read_csv(NFW_SUMMARY_PATH)
    nfw_df.columns = [c.strip() for c in nfw_df.columns]
    nfw_df['_norm'] = nfw_df['Galaxy'].apply(norm_name)
    print(f"  NFW summary entries: {len(nfw_df)}")

    # Merge on normalized name
    merged = pd.merge(
        h2_df[['galaxy', 'regime', 'delta_sigma_kms', '_norm']].rename(
            columns={'galaxy': 'h2_galaxy', 'regime': 'h2_regime',
                     'delta_sigma_kms': 'h2_delta_sigma_dex'}),
        nfw_df[['Galaxy', 'Regime', 'V200_kms', 'C200', 'rs_kpc',
                'sigma_baseline', 'n_inner_points', 'inner_region_defined',
                'max_abs_delta_sigma', 'V_bar_over_V_NFW_inner_median', '_norm']].rename(
            columns={'Galaxy': 'nfw_galaxy', 'Regime': 'nfw_regime'}),
        on='_norm', how='inner'
    )
    print(f"  Matched galaxies: {len(merged)}")

    # Filter to galaxies with valid inner regions
    valid = merged[merged['inner_region_defined'] == True].copy()
    valid = valid.dropna(subset=['max_abs_delta_sigma'])
    print(f"  With valid inner region (n_inner >= 3): {len(valid)}")

    # Classify suppression
    valid['nfw_suppression_class'] = valid['max_abs_delta_sigma'].apply(classify_suppression)
    valid['h2_abs_delta_sigma_dex'] = valid['h2_delta_sigma_dex'].abs()

    # Build output dataframe
    out_cols = ['h2_galaxy', 'h2_regime', 'h2_delta_sigma_dex',
                'V200_kms', 'C200', 'rs_kpc', 'sigma_baseline',
                'n_inner_points', 'max_abs_delta_sigma',
                'V_bar_over_V_NFW_inner_median', 'nfw_suppression_class']
    out_df = valid[out_cols].copy()
    out_df.to_csv(OUT_CSV, index=False, float_format='%.6f')
    print(f"\nSaved: {OUT_CSV}")

    # Spearman correlation: max|Dsigma_NFW| vs Vbar/VNFW_inner_median
    spearman_result = None
    corr_data = valid.dropna(subset=['max_abs_delta_sigma', 'V_bar_over_V_NFW_inner_median'])
    if len(corr_data) >= 5:
        rho, pval = stats.spearmanr(corr_data['V_bar_over_V_NFW_inner_median'],
                                     corr_data['max_abs_delta_sigma'])
        spearman_result = (rho, pval, len(corr_data))

    # By-regime statistics
    regime_stats = {}
    for regime in ['baryon-dom', 'balanced', 'DM-dom']:
        sub = valid[valid['h2_regime'] == regime]
        if len(sub) == 0:
            continue
        ds = sub['max_abs_delta_sigma'].dropna()
        regime_stats[regime] = {
            'N': len(sub),
            'mean': float(ds.mean()),
            'median': float(ds.median()),
            'std': float(ds.std()),
            'min': float(ds.min()),
            'max': float(ds.max()),
        }

    # Suppression classification counts
    suppression_counts = valid['nfw_suppression_class'].value_counts().to_dict()

    # Consistency assessment
    # H2 shows delta_sigma_dex ~ 0 (null result: scatter does not change with perturbation)
    # NFW: if max|Dsigma_NFW| is small, it's consistent with H2 null result
    n_strong   = suppression_counts.get('strong_suppression', 0)
    n_moderate = suppression_counts.get('moderate', 0)
    n_weak     = suppression_counts.get('weak_suppression', 0)
    n_total    = len(valid)

    h2_null_consistent = (n_strong + n_moderate) / n_total if n_total > 0 else 0

    # Overall median of max|Dsigma_NFW|
    overall_median = float(valid['max_abs_delta_sigma'].median())
    overall_mean   = float(valid['max_abs_delta_sigma'].mean())

    # ---- Write narrative report ----
    with open(OUT_REPORT, 'w') as rpt:
        rpt.write("H2 vs NFW Comparison Report\n")
        rpt.write("=" * 60 + "\n\n")

        rpt.write("OVERVIEW\n")
        rpt.write("-" * 40 + "\n")
        rpt.write(f"H2 fleet galaxies analyzed:           {len(h2_df)}\n")
        rpt.write(f"NFW perturbation results available:   {len(nfw_df)}\n")
        rpt.write(f"Matched and valid (n_inner >= 3):     {len(valid)}\n\n")

        rpt.write("H2 NULL RESULT CONTEXT\n")
        rpt.write("-" * 40 + "\n")
        rpt.write("H2 (Adaptive Kernel Framework) finds inner-region delta_sigma ~ 0\n")
        rpt.write("in the operational log-residual RMS metric used by its diagnostic\n")
        rpt.write("implementation, indicating scatter-neutral behavior under bounded\n")
        rpt.write("kernel adaptation. This is the reference pattern for comparison.\n\n")

        rpt.write("NFW PERTURBATION SENSITIVITY RESULTS\n")
        rpt.write("-" * 40 + "\n")
        rpt.write("Perturbations applied: V200 +-10%/+-20%, C200 +-10%/+-20%, both jointly\n")
        rpt.write("Scatter metric: RMS of log10(V_tot) - log10(V_obs) in inner region\n")
        rpt.write(f"Inner region: R < 0.5 * R_max\n\n")

        rpt.write(f"Overall max|Delta_sigma_NFW| [dex]:\n")
        rpt.write(f"  Mean:   {overall_mean:.4f} dex\n")
        rpt.write(f"  Median: {overall_median:.4f} dex\n")
        rpt.write(f"  (N = {n_total})\n\n")

        rpt.write("Suppression classification (max|Delta_sigma|):\n")
        rpt.write(f"  < 0.01 dex (strong suppression):    {n_strong:3d} / {n_total} "
                  f"({100*n_strong/n_total if n_total else 0:.1f}%)\n")
        rpt.write(f"  0.01-0.05 dex (moderate):           {n_moderate:3d} / {n_total} "
                  f"({100*n_moderate/n_total if n_total else 0:.1f}%)\n")
        rpt.write(f"  > 0.05 dex (weak suppression):      {n_weak:3d} / {n_total} "
                  f"({100*n_weak/n_total if n_total else 0:.1f}%)\n\n")

        rpt.write("BY REGIME\n")
        rpt.write("-" * 40 + "\n")
        for regime, stats_d in regime_stats.items():
            rpt.write(f"  {regime} (N={stats_d['N']}):\n")
            rpt.write(f"    mean={stats_d['mean']:.4f}, median={stats_d['median']:.4f}, "
                      f"std={stats_d['std']:.4f}\n")
            rpt.write(f"    range=[{stats_d['min']:.4f}, {stats_d['max']:.4f}]\n")
        rpt.write("\n")

        if spearman_result is not None:
            rho, pval, n_corr = spearman_result
            rpt.write("SPEARMAN CORRELATION\n")
            rpt.write("-" * 40 + "\n")
            rpt.write(f"max|Delta_sigma_NFW| vs V_bar/V_NFW_inner (N={n_corr}):\n")
            rpt.write(f"  Spearman rho = {rho:.3f}, p-value = {pval:.4e}\n")
            if pval < 0.05:
                rpt.write(f"  -> Statistically significant correlation (p < 0.05)\n")
                if rho > 0:
                    rpt.write("  -> Baryon-dominated galaxies show larger NFW parameter sensitivity\n")
                else:
                    rpt.write("  -> DM-dominated galaxies show larger NFW parameter sensitivity\n")
            else:
                rpt.write("  -> No statistically significant correlation (p >= 0.05)\n")
            rpt.write("\n")

        rpt.write("CONSISTENCY WITH H2 NULL RESULT\n")
        rpt.write("-" * 40 + "\n")
        rpt.write(f"H2 reports delta_sigma ~ 0 in the same log-residual scatter metric\n")
        rpt.write(f"for all regime types (scatter-neutral behavior under bounded perturbation).\n")
        rpt.write(f"NFW analysis: {100*h2_null_consistent:.1f}% of galaxies show\n")
        rpt.write(f"  max|Delta_sigma_NFW| < 0.05 dex (strong or moderate suppression).\n\n")

        if overall_median < 0.01:
            assessment = "CONSISTENT: NFW fits show strong scatter suppression (max|Dsigma| << 0.01 dex median), broadly aligned with H2 null result."
        elif overall_median < 0.05:
            assessment = "BROADLY CONSISTENT: NFW fits show moderate scatter sensitivity (0.01-0.05 dex median). NFW perturbations cause somewhat larger scatter changes than H2 kernel perturbations, but still relatively small."
        else:
            assessment = "INCONSISTENT: NFW fits show substantial scatter sensitivity (>0.05 dex median). Within the tested perturbation ranges, NFW parameter variations produce larger inner-region scatter changes than the H2 null-response pattern."

        rpt.write(f"Overall assessment:\n  {assessment}\n\n")

        rpt.write("METRIC COMPARABILITY NOTE\n")
        rpt.write("-" * 40 + "\n")
        rpt.write("This comparison assumes that the H2 operational inner-region scatter\n")
        rpt.write("diagnostic is the same log10-based RMS metric used internally in the\n")
        rpt.write("H2 diagnostic implementation. Under that assumption, both H2 and NFW\n")
        rpt.write("results are being interpreted in the same dex/log-residual metric.\n")
        rpt.write("The legacy fleet-summary column name 'delta_sigma_kms' is therefore\n")
        rpt.write("treated here as a mislabeled field rather than a distinct physical unit.\n")

    print(f"Saved: {OUT_REPORT}")
    print("\nComparison report summary:")
    print(f"  Matched galaxies with valid inner region: {len(valid)}")
    print(f"  Overall median max|Dsigma_NFW|: {overall_median:.4f} dex")
    print(f"  Strong suppression (<0.01 dex): {n_strong}/{n_total}")
    print(f"  Moderate (0.01-0.05 dex): {n_moderate}/{n_total}")
    print(f"  Weak (>0.05 dex): {n_weak}/{n_total}")

    print("\nFirst 10 rows of comparison summary:")
    print(out_df.head(10).to_string())

    return out_df


if __name__ == '__main__':
    main()
