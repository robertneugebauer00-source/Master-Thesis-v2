# LIONESS/BONOBO single-site networks — final commit package

This folder contains everything added to the `sgh-chemical-networks` pipeline
(Master-Thesis-v2) for the single-site (per-site) network extension, laid out
**repo-relative** so it can be copied over the repository root and committed
as-is. It corresponds to branch `lioness-bonobo` (5 commits), also provided
as `lioness_bonobo_branch.patch` for `git am`.

## What this adds

| Step | File | What it does |
|---|---|---|
| — | `R/05_lioness.R` | Single-sample network estimators: `lioness_cor`, `bonobo_cor` (netZooPy-calibrated), `lioness_glasso`, plus `single_site_summary` and caching helpers. |
| S6 | `analysis/S6_single_networks.Rmd` | Builds one network per site for every unit × data kind, summarises each on the aggregate StARS skeleton, writes `master_single_networks.csv`. |
| S7 | `analysis/S7_site_level_tests.Rmd` | Site-level stress-gradient tests: pooled inference (sign count, median ρ, Stouffer) with a within-unit permutation null; circularity-free UDF proxy; trivial-score benchmark; biomarker F-test. |
| S8 | `analysis/S8_backbone_sensitivity.Rmd` | Backbone sensitivity: repeats the S7 tests with `lioness_glasso` (per-site glasso refits at the fixed aggregate λ) so the S7 sign pattern is checked against a different backbone. Requires one S6 knit with `LIONESS_GLASSO <- TRUE`. |
| S9 | `analysis/S9_simulation.Rmd` | Paired chemical + microbial simulation with known ground truth, calibrated on the real Elbe unit: per-site decoupling κ (intact vs broken correlation structure, rising with stress), Dirichlet-multinomial community with amount/structure response channels, recovery of κ by the unmodified S6 estimators, the S7 biomarker F-test replayed under H0/H_amount/H_structure/H_both, and a power grid over n × β2. Designs the future same-site study. |
| — | `tests/lioness_test.R` | 23 unit checks for the estimators (LIONESS mean identity, BONOBO PSD/bounds, netZooPy cross-check to 1e-6, skeleton restriction). |
| — | `data-raw/Danube_TRIDENT_mapping_reconstructed.csv` | Compound→SMILES mapping reconstructed via InChIKey (needed for TU/EC10). |
| — | `run_all.R` | Pipeline runner, now with `S6`, `S7`, `S8`, `S9` steps (all `FALSE` by default). |

## Verification

**Every file on this branch was executed in a real R 4.5.1 environment on the
real data before delivery** — see `VERIFICATION.md`. Short version: test suite
23/23 PASS (incl. live netZooPy cross-check to 1e-15), S9 full run matches the
independent Python mirror, S6/S7/S8 full code-path runs produce all tables and
figures, the 7-commit patch applies cleanly to a fresh clone of the base.
Three bugs found by execution were fixed (commit `ba876ce`), and one
pre-existing mainline fragility (StARS constant-column aborts on units with
rare compounds) is documented there.

## How to commit

```bash
# option A: copy files over the repo root, review, commit
cp -r final_commit/* /path/to/Master-Thesis-v2/   # results/ optional, see below
git add -A && git commit

# option B: apply the whole branch as patches
cd /path/to/Master-Thesis-v2 && git am < lioness_bonobo_branch.patch
```

`results/` holds the analysis outputs (tables + figures) for reference; they
are reproducible by knitting S6–S8 and need not be committed. `results/s8`
additionally contains the exact Python replication of the S8 run
(`prep.py`, `s8_run.py`, `s8_stats.py`, logs) used to produce the delivered
numbers, in case the R run and the delivered CSVs should be compared.

## How to run

```r
setwd("repo root"); source("run_all.R")     # S1–S8 off by default; flip STEPS
# S8 needs S6 knitted once with LIONESS_GLASSO <- TRUE (costs n glasso
# solves per unit; ~17 units x ~50 sites, minutes per unit)
```

## Headline results (delivered numbers)

- **S6**: 34/34 unit×kind combos produced per-site networks; BONOBO per-unit
  Spearman ρ(strength, sum-TU) negative in 13/17 units (median −0.15).
- **S7**: pooled Stouffer Z = −2.56 (permutation p = 0.005) for BONOBO vs
  Z = +3.76 for LIONESS — the estimators bracket and disagree; the external
  UDF proxy sides with BONOBO (negative in 3/3 basins, Elbe p = 0.02);
  Elbe survives trivial-score controls (partial ρ = −0.303); biomarker test
  NEGATIVE (strength adds nothing beyond sum-TU for UDF, p = 0.98).
- **S8** (backbone sensitivity): the sign of the site-level coupling is
  backbone-dependent. On the glasso backbone, lioness_glasso is strongly
  POSITIVE (Stouffer Z = +15.5; Elbe UDF ρ = +0.39, p < 1e-4) — but its
  per-basin partial ρ given the trivial scores collapses to ~0 or negative
  in 4/5 basins, i.e. the glasso-backbone strength is largely a
  detection-richness score. BONOBO recomputed on the same skeleton keeps its
  negative sign (Z = −1.66, attenuated vs S7's −2.56 — direction robust,
  magnitude skeleton-sensitive). Per-site agreement between the two
  backbones is essentially zero (Spearman 0.09). Net effect on the thesis:
  the correlation-backbone BONOBO result is the only site-level signal that
  survives trivial-score controls; report the coupling as hypothesis-
  generating and estimator- AND backbone-dependent. See
  `results/s8/tableS8_vs_S7_comparison.csv` and the S8 Rmd's reading guide.
- **S9** (simulation, validated by Python mirror in `results/s9/`): on known
  ground truth, BONOBO strength recovers the true decoupling κ
  (Spearman ρ = −0.47…−0.50, AUC 0.73–0.76 across n = 50–500) while LIONESS
  points the wrong way (ρ ≈ +0.11, AUC < 0.5 — the leverage conflation
  reproduced in silico). The S7 biomarker F-test is well calibrated
  (H0 false-positive rate 0.0–0.1), specific (H_amount rejection 0.10) and
  sensitive (H_both 0.80 at n = 200). Power for a structure-driven biological
  response: 0.3 at n = 50 → 0.9 at n = 200 → 1.0 at n = 500 for β2 = 0.3;
  ≥ 0.9 already at n = 200 for β2 ≥ 0.6. Read: a future same-site
  chemical + eDNA study should aim for n ≳ 200 sites. Caveat: identifiability
  study — it shows the pipeline CAN find structure-driven biology at these
  sample sizes, not that such biology exists.

## Setting up on a new machine

1. Install R 4.5.x and RStudio (the lockfile was made with R 4.5.3).
2. Clone or copy this folder anywhere -- no paths need editing; the scripts locate the repo themselves.
3. Open `sgh-chemical-networks.Rproj`. renv bootstraps itself via `.Rprofile`.
4. Run `renv::restore()` once to install the exact package versions from `renv.lock`
   (CRAN, Bioconductor and the GitHub packages NetCoMi, SpiecEasi, SPRING).
5. `source("run_all.R")`. Inputs are read from `data-raw/`, results go to `outputs/`.

After installing or updating a package, run `renv::snapshot()` and commit `renv.lock`.
