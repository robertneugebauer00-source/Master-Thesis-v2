# Verification report — lioness-bonobo branch (7 commits)

Everything below was executed in a real R environment (R 4.5.1, Posit build,
extracted user-locally; packages compiled from source: igraph 2.3.3, huge 2.0.1,
WGCNA 1.74, tidyverse, readxl, uwot, knitr, …) on the real uploaded data
(Finckh_Carmona_2023_PANGAEA_R1.xlsx, MDL_CECs_summary.csv,
TRIDENT_prediction_results_EC10_GRO.csv, reconstructed Danube mapping).
Date: 17–18.09.2026.

## What was executed, and what happened

| Component | How it was verified | Result |
|---|---|---|
| `tests/lioness_test.R` | Full run in R | **23/23 PASS, 0 skips** — incl. test 5 (`lioness_glasso` with huge) and the live Python-netZooPy cross-check (max diff 1.4e-15) |
| Foundation (`00_setup`, `01_data`, `02_functions`, `05_lioness`) | Sourced with real data | Units, node sets and sum-TU match the independently validated Python prep exactly (e.g. EUS_029 sum-TU = 275.7981; Elbe 153 sites × 182 compounds; MDL match 608/610) |
| `S9_simulation.Rmd` | **Full run** (purl mode), all chunks | All tables + both figures written. Results match the independent Python mirror within sampling noise (see below) |
| `S6_single_networks.Rmd` | Functional run on 2 basins (Elbe + Danube), all 3 estimators incl. `LIONESS_GLASSO <- TRUE` | Caches, master table (840 rows), per-unit figures all produced; per-unit `tryCatch` resilience demonstrated |
| `S7_site_level_tests.Rmd` | Full code-path run (NPERM = 200 in the test copy) | All 5 tables + 4-panel figure produced; pooled models, biomarker test, trivial-score benchmark all execute |
| `S8_backbone_sensitivity.Rmd` | Full code-path run (NPERM = 200) | All 5 tables + 2-panel comparison figure produced |
| All Rmds + R files + `run_all.R` | `parse()` on every file / purl'd chunk | All parse clean |

## Bugs found by execution and fixed (commit `ba876ce`)

1. **S9 — NA in the stress pool.** 11 Elbe sites have no TU column (no EC10
   coverage), so `sum_tu[elbe_sites]` contained NAs and the resampled stress
   axis poisoned the simulator. Fixed: `TU_POOL <- TU_POOL[!is.na(TU_POOL)]`.
2. **S7/S8 — single-unit `replicate()` simplification.** With one unit, the
   permutation null collapsed to a vector and `rowMeans` failed. Fixed by
   forcing the matrix shape (a no-op for the full 17-unit run).
3. **S7 — hardcoded `seq_len(17)`** in the forest plot. Now
   `seq_along(units17)` so subset runs don't crash the figure.

## S9: R run vs independent Python mirror

| Metric | R (this run) | Python mirror |
|---|---|---|
| ρ(BONOBO strength, true κ) | −0.45 … −0.49 | −0.47 … −0.50 |
| ρ(LIONESS strength, true κ) | +0.13 … +0.14 | +0.09 … +0.15 |
| AUC BONOBO / LIONESS | 0.71–0.76 / 0.42–0.43 | 0.73–0.76 / 0.42–0.45 |
| H0 false-positive rate | 0.00–0.10 | 0.00–0.10 |
| Power, β2 = 0.3 (n = 50/100/200/500) | 0.30 / 0.80 / 0.90 / 1.00 | 0.30 / 0.60 / 0.90 / 1.00 |
| H_amount control (specificity) | 0.00 | 0.10 |
| H_both control (sensitivity) | 0.85 | 0.80 |

## Test-environment deviations (deliverable code unaffected)

- `magick` stub package: `library(magick)` appears in `00_setup.R` but magick
  is **never called** anywhere in `R/` or `analysis/` (verified by grep); the
  system lacks ImageMagick, so an empty stub satisfied the loader.
- S6 smoke used `MIN_DETECT_UNIT <- 5L` and 2 basins; S7/S8 test copies used
  NPERM = 200. Full 17-unit S6–S8 runs on the user's machine use the canonical
  caches and parameters.
- `pagedown` not installed (optional, guarded in `00_setup.R`).

## Pre-existing issue discovered (NOT introduced by this branch)

`fit_glasso()` / huge-StARS can fail with *"Raw data x contains a constant
column"* on units whose node set contains compounds detected at very few sites
(MIN_DETECT_UNIT = 2): a StARS subsample that misses all detections of such a
compound makes it constant, and huge aborts. Observed on Basin_Danube,
Basin_Rhine, Basin_Tajo, Basin_Basque in this environment (Elbe passed). This
is seed- and huge-version-dependent and can affect re-knits of
`03_mainline.Rmd` on fresh machines. The per-unit `tryCatch` contains it, but
affected units silently lose their network. Recommended (not applied, mainline
behaviour change): raise `MIN_DETECT_UNIT` or drop within-subsample constant
columns in `fit_glasso`.

## How to push

```bash
cd /path/to/Master-Thesis-v2
git fetch origin
git push origin lioness-bonobo
# or, from a fresh clone, apply the whole branch:
git am < lioness_bonobo_branch.patch   # 7 commits, 19cadf9..ba876ce
```
