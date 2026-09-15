# sgh-chemical-networks

Chemical co-occurrence networks in European rivers, used to test the **Stress Gradient
Hypothesis** (SGH): does modular network organisation decline along a stress gradient?

Measured concentrations of ~600 compounds at sites in the Elbe, Rhine and Danube are
turned into a sparse **partial-correlation network** per spatial unit (graphical lasso,
StARS-selected), partitioned into modules, and tested against degree-preserving null
models. Toxicological risk is overlaid as a second, independent node attribute.

MSc thesis, Robert Neugebauer. Supervisor: Pedro Inostroza. Builds on P. A. Inostroza's
vPI3 workflow.

## Layout

```
run_all.R                          run everything, in order, with timing and a log
_targets.R                         OPTIONAL {targets} pipeline (caching/parallelism)
R/
  00_setup.R                       palette, packages, paths, ALL analysis parameters
  01_data.R                        loading, MDL/2 floor, toxic units, spatial units
                                   (country groups, basins, sections, reaches)
  02_functions.R                   prep / fit / metrics / modules / risk / plotting
                                   + fitted-artefact cache, run manifest,
                                   module-stability bootstrap, signed-edge check
  null_models_vRN.R                the null-model hierarchy (ER, config, Chung-Lu,
                                   generative families, and the data-level pipeline null)
  highlights_vRN.R                 tag index: hl(), hl_jump(), hl_sections()
  99_utils.R                       figure refresh + standalone HTML report builder
                                   (side-effect-free on source; call main())
analysis/
  03_mainline.Rmd                  THE PIPELINE -- produces Tables 1-4 and Figures 1-4,
                                   and saves fitted_<kind>_<unit>.rds per unit
  S1_method_justification.Rmd      prevalence sweep, node-set audit, glasso vs MB, MDL/2,
                                   unit inventory, small-unit StARS depth, signed edges
  S2_alternative_communities.Rmd   SBM / DC-SBM, ICL-selected, no random-graph null
  S3_crossbasin.Rmd                rarefaction to common n, module sharing, basin panel
  S4_null_hierarchy.Rmd            all ten null families + the data-level pipeline null
  S5_stress_gradient.Rmd           the SGH test itself: modularity along a stress index
tests/smoke_test.R                 end-to-end test on synthetic data (no private data)
data-raw/fetch_pangaea.R           download the PANGAEA dataset by DOI
archive/legacy_null_test.R         superseded ER-only null, kept for reproducibility
docs/split_map.md                  how this was split out of the monolithic script
data/README.md                     what data is needed and where to point the code
```

## What changed on 15.09.2026 (consolidation pass)

- `run_unit()` now **saves the fitted networks** (`fitted_<kind>_<unit>.rds`: adjacency
  matrices, memberships, StARS lambdas). The figure refresh in 99_utils.R and the S2/S3
  companions read these instead of refitting -- "figure refresh without refitting" is now
  literally true, and S2 no longer needs the external `_fit_cache.rds` from
  `vRN_followups.R` (a legacy fallback is kept).
- **Single sources of truth**: the metric-row builder (`metric_row()`) lives once in
  02_functions.R (the two old copies had already drifted); reach units are defined once
  in 01_data.R (was: 99_utils.R *and* an external script).
- `null_test_all()` is deprecated (it was never called); the generative null families and
  the pipeline null are now wired into **S4** instead of sitting dormant.
- **S5** adds the explicit SGH test: a per-unit stress index from toxic units vs
  modularity, with n as covariate and a module-stability bootstrap.
- Each run writes `run_manifest.json` (parameters, package versions, git commit) into the
  output tree, so every number is traceable to its configuration.
- 99_utils.R is side-effect-free when sourced (entry point: `main()`; run_all.R drives it).
- Repo-root detection no longer depends on a stray `sp` variable.
- `tests/smoke_test.R` + GitHub Actions run the core pipeline on synthetic planted-module
  data -- no private data needed.

The split rule: **mainline is whatever produces Tables 1-4 and Figures 1-4;
supplementary is whatever defends a choice.**

## Quick start

```r
Sys.setenv(SGH_PROJ = "C:/path/to/your/R")   # where the data lives
setwd("path/to/sgh-chemical-networks")
source("run_all.R")                          # edit STEPS at the top to pick what runs
```

Or load the foundation and work interactively:

```r
source("R/00_setup.R"); source("R/01_data.R"); source("R/02_functions.R")
```

Verify the checkout runs (synthetic data, no downloads, needs only igraph + huge):

```r
source("tests/smoke_test.R")   # or: Rscript tests/smoke_test.R
```

## Reproducibility

Package versions are not yet pinned by a lockfile. To create one after installing
the dependencies, run once from the repo root:

```r
install.packages("renv"); renv::init()   # writes renv.lock -- commit it
```

Every pipeline run already writes `run_manifest.json` (parameters, package
versions, git commit) into the output tree.

## Method in one paragraph

Two chemicals correlate across sites mostly because some sites are polluted and some
are clean. A **partial** correlation asks whether they still co-vary once every other
compound is accounted for, so the shared pollution gradient does not manufacture edges.
With more compounds than sites the partial-correlation matrix cannot be inverted
directly, so the **graphical lasso** estimates a regularised version and **StARS**
picks the penalty by edge stability across site resamples. Modules are found by greedy
modularity maximisation; an edge-level test over ~18,500 compound pairs per basin shows
they are **emission-source** compartments, not mechanistic ones.

## On the null model

Newman's Q already contains a configuration-model expectation (`k_i*k_j/2m`), so an
external null ensemble must preserve degrees to test the same hypothesis the statistic
assumes. `null_test_multi()` therefore reports **config** and **Chung-Lu** first and
Erdos-Renyi last. This is not cosmetic: **8 of 17 units change verdict depending on the
null, and 5 flip sign.** Quote config in the main text; keep ER for average path length,
where it is the conventional reference.

Change the null in one place -- the `nulls =` argument in `analysis/03_mainline.Rmd`.

## Finding your way around the code

```r
source("R/highlights_vRN.R")
hl()                 # index of @KEY / @NULL / @PARAM / @DECISION / @CAVEAT markers
hl_jump("NULL", 2)   # jump straight to the null test
```

## Key parameters

All in `R/00_setup.R`:

| Parameter | Value | Controls |
|---|---|---|
| `PREV_THRESH` | 0.10 | prevalence filter -- which compounds become nodes |
| `NODESET_SCOPE` | `global` | filter applied once over all reference sites |
| `STARS_REPNUM` | 100 | StARS subsampling depth (measured, not guessed) |
| `RQ_THRESH` | 0.02 | absolute-risk-driver threshold |
| `BASIN_MIN_N` / `SECTION_MIN_N` | 20 / 30 | minimum sites per unit |

## Data

Not included -- see [`data/README.md`](data/README.md).

## Licence

Code is MIT licensed (see `LICENSE`). The datasets are **not** covered by it and remain
under their original terms.
