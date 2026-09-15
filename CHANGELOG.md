# Changelog

## 15.09.2026 — consolidation + expansion pass

Applied the recommendations from the repository review (`sgh_repo_review_and_roadmap.md`).
**No analysis parameter, seed, estimator or filter value was changed** — every pre-existing
number the pipeline produces is untouched. Changes are plumbing, deduplication, and new
(optional) analyses.

### Fixed / consolidated

- **`metric_row()`** — the metric-row builder was defined twice (`03_mainline.Rmd` and
  `99_utils.R`) and the copies had drifted (the 99_utils copy lacked `Nodes (pool)`).
  Now a single definition in `R/02_functions.R`; both call sites use it.
- **Reach units** — were defined in `99_utils.R` *and* in the external `vRN_followups.R`.
  The single definition now lives in `R/01_data.R` (`make_reach()` / `reach_units`).
- **`99_utils.R` is side-effect-free** — sourcing it defines functions only; call
  `refresh_all_figures()` / `main()` (run_all.R's `utils` step does this for you; a
  `sys.nframe()` guard keeps `Rscript R/99_utils.R` working).
- **Refresh without refitting, for real** — `refresh_unit_figures()` reuses the cached
  `fitted_*.rds` when present and only refits units without a cache (or with a
  glasso-only cache from the optional targets pipeline — detected and refit).
- **Repo-root resolution** — `00_setup.R` no longer trusts a pre-existing `sp` variable;
  it validates it, else walks up to `run_all.R`. `highlights_vRN.R` paths are root-aware
  and `hl_jump()` degrades gracefully outside RStudio.
- **`null_test_all()` deprecated** — it was dead code; `.Deprecated()` points to
  `null_test_multi()`.
- **S2 self-contained** — reads the mainline's fitted rds instead of requiring the
  external `_fit_cache.rds` from `vRN_followups.R`. The legacy cache is merged in as a
  base layer so units the mainline never runs (the four `Reach_*` units) keep their
  coverage. Caveat: if your legacy cache was fitted with different settings, re-knitting
  S2 changes its numbers — the mainline rds are canonical.
- **S3 basin panel** reads cached basin fits instead of refitting
  (`BASIN_CMP_REFIT <- TRUE` still forces a refit).
- Docs updated: `docs/split_map.md` (run instructions for the repo layout, to-do items
  closed), `README.md`, `data/README.md`.

### Added

- **Fitted-artefact cache** — `run_unit()` saves `fitted_<kind>_<unit>.rds` per unit
  (adjacency matrices, memberships, StARS lambdas). Helpers: `kind_dir()`,
  `fitted_path()`, `load_fitted()`, `save_fitted()`.
- **Run manifest** — `write_run_manifest()` writes `run_manifest.json` (parameters,
  package versions, git commit) into the output tree every mainline run.
- **`analysis/S4_null_hierarchy.Rmd`** — wires in the previously dormant null machinery:
  `null_family_test()` + `plot_null_grid()` for all ten generative families, and the
  data-level `pipeline_null(type = "grad")` with density-matched comparison, for the
  three whole-basin networks.
- **`analysis/S5_stress_gradient.Rmd`** — the explicit SGH test: per-unit stress index
  from toxic units vs glasso modularity (Spearman + weighted lm, with n as covariate),
  a within-basin Upper→Lower gradient panel, and a module-stability bootstrap
  (`bootstrap_module_stability()`, which compares against the cached mainline partition).
- **S1 additions** — automated unit inventory across prevalence thresholds (the open
  `PREV_THRESH = 0.05` question), a StARS-depth sweep for the small section units, and
  the signed-edge check (`signed_edge_check()`).
- **`tests/smoke_test.R` + `.github/workflows/smoke.yaml`** — end-to-end test on
  synthetic planted-module data; needs only igraph + huge, no private data.
- **`_targets.R`** — optional {targets} pipeline (caching + unit-level parallelism).
  `run_all.R` remains the canonical entry point.
- **`data-raw/fetch_pangaea.R`** — downloads the concentration dataset by DOI
  (10.1594/PANGAEA.960272).

### Deliberately not done

- `renv.lock` — needs an R environment with the packages installed; run
  `renv::init()` once on your machine (see README → Reproducibility).
- ClassyFire cache — requires live API calls; the S1 chunk stays `eval = FALSE` by design.
- Full R-package conversion (DESCRIPTION/R CMD check) — deferred; the smoke test +
  CI give the same safety net without restructuring the repo.
