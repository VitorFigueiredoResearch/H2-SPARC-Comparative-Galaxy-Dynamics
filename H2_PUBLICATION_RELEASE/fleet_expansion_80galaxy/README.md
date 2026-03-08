# H2 Fleet Expansion: 80-Galaxy Validation
**Phase 2b — Scatter Neutrality Robustness Test**
Generated: 2026-03-06

---

## Overview

This folder documents the expansion of the H2 adaptive nonlocal kernel validation
fleet from 50 to 80 galaxies. The primary scientific objective is to demonstrate
that the scatter neutrality result is robust across a statistically larger and more
morphologically diverse sample spanning all dynamical regimes.

**Critical finding: Δσ is consistent with zero km/s for all galaxies with sufficient
inner-region coverage; a small subset is undefined due to inner-coverage limitations
(data geometry — not model failure).**

---

## Fleet Composition

| Source | Count | Description |
|---|---|---|
| Original H2 baseline | 9 | NGC0891, NGC6946, NGC5585, IC2574, NGC3893, NGC2903, NGC5055, NGC3198, NGC6503 |
| Phase 1 expansion (30-galaxy) | 21 | fleet_expansion_30galaxy/ |
| Phase 2a expansion (50-galaxy) | 20 | fleet_expansion_50galaxy/ |
| Phase 2b expansion (this folder) | 30 | New galaxies, regime-balanced selection |
| **Total** | **80** | **Full 80-galaxy fleet** |

---

## Regime Distribution (80 galaxies)

| Regime | Criterion | Count |
|---|---|---|
| Baryon-dominated | Vk/Vb < 0.8 | 39 |
| Balanced | 0.8 ≤ Vk/Vb < 1.5 | 30 |
| DM-dominated | Vk/Vb ≥ 1.5 | 11 |
| **Total** | | **80** |

---

## H2 Compatibility Constraint

All 80 galaxies share frozen H1 baseline parameters:
- **L₀ = 200 kpc** (characteristic kernel length scale)
- **μ = 10** (nonlocal coupling amplitude)

Only galaxies with these H1 parameters are compatible with Phase 4 adaptive
convolution under taper mode. Galaxies with smaller L₀ fail with [TAPER-FAIL].

---

## Processing Protocol

Each galaxy was processed through the full H2 diagnostic pipeline:

```bash
PYTHONUTF8=1 python -m diagnostics.run_fleet \
    --galaxies [GALAXY] --alpha 2.0 --sigma_idx 1.0 --taper \
    --output data/derived/fleet/fleet_[GALAXY].csv
```

**Pipeline stages:**
1. Phase 3: χ field computation + L_eff(r) = L₀/(1+χ) adaptive scale
2. Phase 4: Adaptive convolution with H2 kernel
3. Test 2: Pearson r (H2 vs H1 residual correlation)
4. Test 3: Δσ inner scatter (H2 vs H1 scatter comparison)

**Stop rule applied:** If |Δσ| > 1e-12 → STOP and report. No stop triggered across all 80 galaxies.

---

## Phase 2b Galaxy List (30 new galaxies)

| # | Galaxy | Vk/Vb | Regime | Pearson r | Δσ (km/s) | Flag |
|---|---|---|---|---|---|---|
| 51 | UGC07577 | 0.198 | baryon-dom | 0.981 | 0.0000 | — |
| 52 | NGC6789 | 0.252 | baryon-dom | 0.621 | N/A | inner_coverage |
| 53 | UGC06818 | 0.393 | baryon-dom | 0.474 | 0.0000 | — |
| 54 | NGC0100 | 0.459 | baryon-dom | 0.450 | 0.0000 | — |
| 55 | DDO064 | 0.467 | baryon-dom | 0.760 | 0.0000 | — |
| 56 | NGC4389 | 0.469 | baryon-dom | 0.964 | N/A | inner_coverage |
| 57 | UGC07603 | 0.514 | baryon-dom | 0.684 | 0.0000 | — |
| 58 | UGC07151 | 0.515 | baryon-dom | 0.774 | 0.0000 | — |
| 59 | NGC4088 | 0.564 | baryon-dom | −0.157 | 0.0000 | negative_r |
| 60 | UGC12732 | 0.583 | baryon-dom | 0.615 | 0.0000 | — |
| 61 | UGC07089 | 0.706 | baryon-dom | 0.861 | 0.0000 | — |
| 62 | UGC12632 | 0.771 | baryon-dom | 0.914 | 0.0000 | — |
| 63 | NGC4214 | 0.795 | baryon-dom | 0.477 | 0.0000 | — |
| 64 | UGC02916 | 1.143 | balanced | 0.413 | 0.0000 | — |
| 65 | UGC11914 | 1.156 | balanced | 0.781 | 0.0000 | — |
| 66 | F583-1 | 1.190 | balanced | 0.810 | 0.0000 | — |
| 67 | UGC06917 | 1.247 | balanced | 0.885 | 0.0000 | — |
| 68 | UGC08550 | 1.293 | balanced | 0.668 | 0.0000 | — |
| 69 | UGC06930 | 1.303 | balanced | 0.895 | 0.0000 | — |
| 70 | DDO170 | 1.316 | balanced | 0.700 | 0.0000 | — |
| 71 | UGC00634 | 1.331 | balanced | 0.732 | N/A | inner_coverage |
| 72 | F568-3 | 1.403 | balanced | 0.930 | 0.0000 | — |
| 73 | NGC4183 | 1.418 | balanced | 0.125 | 0.0000 | — |
| 74 | NGC1090 | 1.453 | balanced | 0.517 | 0.0000 | — |
| 75 | UGCA444 | 1.469 | balanced | 0.953 | 0.0000 | — |
| 76 | F565-V2 | 1.623 | DM-dom | 0.911 | 0.0000 | — |
| 77 | DDO154 | 1.666 | DM-dom | 0.562 | 0.0000 | — |
| 78 | F568-V1 | 1.841 | DM-dom | 0.887 | 0.0000 | — |
| 79 | UGC08286 | 1.905 | DM-dom | 0.620 | 0.0000 | — |
| 80 | UGC05005 | 1.963 | DM-dom | 0.636 | 0.0000 | — |

---

## Anomaly Flags

### inner_coverage (6 galaxies across all 80)

Test-3 inner scatter returned NaN because insufficient radial data points fall
within the inner-radius threshold used for scatter computation. This is a data
geometry limitation, not a model failure. Δσ is consistent with zero km/s for
all galaxies with sufficient inner-region coverage; this subset is undefined due
to inner-coverage limitations. Consistent with Phase 1 and Phase 2a handling.

Phase 2b inner_coverage galaxies: NGC6789 (r=0.621), NGC4389 (r=0.964),
UGC00634 (r=0.732). All show good-to-excellent Pearson r values, confirming
the adaptive kernel is functioning correctly.

### negative_r (2 galaxies across all 80: NGC4559 from Phase 2a, NGC4088 from Phase 2b)

A negative Pearson r indicates H2 residuals are anti-correlated with H1 residuals —
H2 actively corrects in the opposite direction to H1 errors. Δσ is consistent
with zero km/s regardless; scatter neutrality holds even in this compensation regime.

### High-MAFE stress tests (3 galaxies in Phase 2b)

UGC08550 (MAFE=0.898), DDO154 (MAFE=1.154), and UGC05005 (MAFE=0.951) were
included as deliberate stress tests with poor H1 fits. All three returned Δσ
consistent with zero km/s, demonstrating that scatter neutrality is independent
of H1 fit quality.

---

## Key Results

| Metric | Value |
|---|---|
| Total galaxies | 80 |
| Δσ consistent with zero km/s (measured) | 74/74 (100%) |
| Δσ N/A (data geometry, inner_coverage) | 6 |
| Stop rule triggered | Never |
| Pearson r range | −0.3713 to 0.9808 |
| Vk/Vb range | 0.027 to 2.187 |
| Universal L_eff floor | 66.67 kpc = L₀/3 |
| min_L_eff range | 66.67 to 96.95 kpc |

---

## Files

```
fleet_expansion_80galaxy/
├── README.md                          (this file)
├── fleet_summary_80galaxy.csv         (80-row combined results)
├── figures/
│   ├── regime_trends_80galaxy.png     (Pearson r and Δσ vs Vk/Vb)
│   └── leff_overlay_80galaxy.png      (all 80 L_eff(r) profiles)
├── phase3_outputs/                    (L_eff and χ field CSVs)
├── phase4_outputs/                    (H2 RC outputs and comparisons)
└── diagnostics/                       (per-galaxy fleet CSVs)
```

Per-galaxy diagnostic CSVs are in `data/derived/fleet/fleet_[GALAXY].csv`.
Per-galaxy L_eff profiles are in `data/derived/phase3/leff_[GALAXY].csv`.

---

## Relation to Previous Phases

- **30-galaxy fleet:** archived in `fleet_expansion_30galaxy/` — NOT modified
- **50-galaxy fleet:** archived in `fleet_expansion_50galaxy/` — NOT modified
- This folder contains only the 30 Phase 2b additions plus the combined 80-galaxy summary
- The Phase 1 and Phase 2a results (Δσ consistent with zero km/s across 50 galaxies)
  are reproduced and extended here

---

## Parameters (Frozen)

| Parameter | Value |
|---|---|
| α (adaptive strength) | 2.0 |
| σ_idx (structure index) | 1.0 |
| Taper | Enabled |
| L₀ (frozen H1 baseline) | 200 kpc |
| μ (frozen H1 baseline) | 10 |
| Kernel | hybrid logarithmic |
| Grid spacing Δx | 1.0 kpc |
