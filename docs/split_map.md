# Clean_Workflow_vRN -- split into mainline and supplementary

Created 14.09.2026 from `../Clean_Workflow_vRN-14.08.26.Rmd` (2755 lines).
The original is untouched. All code here is **verbatim** -- verified line by line:
every one of the 2171 code lines in the original appears exactly once in this folder,
and the only lines added are the file headers and the `source()` block.

> **Note (15.09.2026):** the verbatim guarantee describes the 14.09.2026 split. The
> consolidation pass of 15.09.2026 changed call sites (not analysis logic) -- see
> "Still to do -- RESOLVED" below for the exact list.

## The rule used for the split

**Mainline** = whatever produces Tables 1-4 and Figures 1-4.  
**Supplementary** = whatever defends a choice.

## Files

| File | Contents | Original chunks |
|---|---|---|
| `00_setup.R` | palette, packages, file paths, **all analysis parameters**, sources `null_models_vRN.R` | palette, package/file setup |
| `01_data.R` | PANGAEA loading, MDL/2 floor, toxic units, spatial unit definitions | data loading, mdl-lookup, TU calc, groups |
| `02_functions.R` | `prep_unit`, `fit_glasso`/`fit_mb`/`fit_wgcna`, metrics + modules, risk drivers, plotting | fn-prep, fn-glasso-mb, fn-wgcna, metric tables, fn-riskdriver, plotting |
| `03_mainline.Rmd` | **the thesis pipeline**: `run_unit()`, run all units, assembly, sanity check | summary, run-all, analysis-pipeline, sanity-check |
| `S1_method_justification.Rmd` | prevalence sweep, node-set audit, netCompare (glasso vs MB), MDL/2 sensitivity, ClassyFire | prev-diagnostic, nodeset-audit, netcompare, mdl-sensitivity, classyfire |
| `S2_alternative_communities.Rmd` | SBM / DC-SBM, ICL-selected, no random-graph null | sbm-blocks |
| `S3_crossbasin.Rmd` | rarefaction to common n, module sharing, Elbe/Danube/Rhine panel | crossbasin, basin fits/panel/barplot |
| `99_utils.R` | figure refresh without refitting, standalone HTML report builder | refresh-figures, build-report |
| `archive/legacy_null_test.R` | the dead ER-only `null_test()` | -- |

## How to run

Updated 15.09.2026 for the repository layout (this section used to show the old
flat OneDrive split folder). From the repo root:

```r
Sys.setenv(SGH_PROJ = "C:/path/to/your/R")   # where the data lives
setwd("path/to/sgh-chemical-networks")
source("run_all.R")                          # edit STEPS at the top to pick what runs
```

or interactively:

```r
source("R/00_setup.R"); source("R/01_data.R"); source("R/02_functions.R")
# then knit analysis/03_mainline.Rmd, or any S1-S5 companion (each sources the three above)
```

## What changed relative to the original

1. **`null_test()` removed** -> `archive/legacy_null_test.R`. It was never called;
   `null_test_multi()` superseded it on 27.08.2026. 56 lines.
2. **`rm(list = ls())` disabled** in `00_setup.R` (line 13). Harmless at the top of a
   knit, destructive when the file is sourced into a live session.
3. **Chunk labels normalised** in the Rmd children: `#RUN ALL NETWORKS; results='hide'`
   was not a valid label and broke `knitr::purl()`. Now `run-all-networks, results='hide'`.
4. Nothing else. No logic, no parameter, no number was touched.

## Still to do -- RESOLVED 15.09.2026

Both items from the original split are now closed:

- `row <- function(m)` was defined twice -- in `run_unit()` (`03_mainline.Rmd`) and in
  `refresh_unit_figures()` (`99_utils.R`) -- and the copies had in fact already drifted
  (the 99_utils copy lacked `Nodes (pool)`). Now a single `metric_row()` in
  `02_functions.R`; both call sites use it.
- The reach units were defined in `99_utils.R` *and* in `vRN_followups.R`. The definition
  now lives in `01_data.R` (`make_reach()` / `reach_units`), next to the other unit
  definitions; the copy in 99_utils.R was removed. `vRN_followups.R` (outside this repo)
  should source it from here.

Further changes on 15.09.2026 (beyond the verbatim split, so the line-for-line
equivalence statement above applies to the 14.09.2026 state):

- `run_unit()` now saves `fitted_<kind>_<unit>.rds` per unit (adjacency matrices,
  memberships, lambdas); 99_utils / S2 / S3 read these instead of refitting. S2 no
  longer needs the external `_fit_cache.rds` from vRN_followups.R (legacy fallback kept).
- `99_utils.R` is side-effect-free on source (`refresh_all_figures()` / `main()`);
  `run_all.R` drives it. The refresh now genuinely skips refitting when cached fits exist.
- `null_test_all()` deprecated (dead code); the generative null families and the pipeline
  null are wired into the new `analysis/S4_null_hierarchy.Rmd`; the explicit SGH test is
  the new `analysis/S5_stress_gradient.Rmd`.
- `00_setup.R` resolves the repo root by walking up to `run_all.R` instead of trusting
  any pre-existing `sp` variable; `highlights_vRN.R` paths are root-aware and `hl_jump()`
  degrades gracefully outside RStudio.

## Tag index

`source("../R/highlights_vRN.R")` then `hl()` -- the `## @KEY / @NULL / @PARAM /
@DECISION / @CAVEAT` markers still work here; point `HL_FILE` at whichever file you want.
