## =============================================================================
## _targets.R -- OPTIONAL {targets} pipeline mirroring run_all.R  [repo 15.09.2026]
##
## STATUS: OPTIONAL. run_all.R remains the canonical entry point -- this file
## exists only for caching (skip refits whose inputs have not changed) and,
## optionally, parallelism. Nothing here changes any analysis number: the same
## prep_unit() / fit_glasso() / greedy_modules() / metrics_from_adj() calls run
## with the same defaults as analysis/03_mainline.Rmd.
##
## HONEST CAVEAT -- this is a BRIDGE, not idiomatic {targets}. The analysis
## code lives in scripts (R/00_setup.R, 01_data.R, 02_functions.R) that create
## globals as a side effect of source(), so the `foundation` target sources
## them and hands the unit grid downstream. Sourced globals do NOT travel
## between targets (each branch runs in a fresh environment), so the unit_fits
## branches re-source the same three scripts -- that costs one data load per
## unit. A full port would wrap prep_unit()/fit_glasso()/greedy_modules() as
## exported functions over explicit data targets (conc_mat, sites_meta, ...)
## and delete the bridge. Until then, prefer correctness over elegance.
##
## RUN:  targets::tar_make()         # from the repo root
##       targets::tar_make(manifest) # just the run manifest
## =============================================================================

library(targets)

tar_option_set(packages = c("igraph", "huge"))

## [repo 15.09.2026] OPTIONAL parallel workers -- commented out on purpose.
## The unit x kind grid (and, in the full pipeline, the N null draws inside
## null_test_multi()) is embarrassingly parallel: units share nothing except
## the read-only data globals. Parallelism is at the UNIT level only, so no
## RNG behaviour inside fit_glasso() / null_test_multi() changes -- each
## branch still sets its own fixed seeds, exactly as in the serial run, and
## results are bit-identical to run_all.R. To enable:
##
##   library(crew)
##   tar_option_set(
##     controller = crew::crew_controller_local(workers = 4)
##   )
##
## (Requires the {crew} package; 4 workers is a sane default for a laptop --
## each branch holds one concentration matrix plus one StARS fit in memory.)

## ---- track the private input file -------------------------------------------
## The pipeline's single irreplaceable input. format = "file" hashes it, so
## downstream targets refit only when the xlsx itself changes. Path logic
## mirrors R/00_setup.R (SGH_PROJ env var, author path as fallback).
tar_target(
  data_files,
  file.path(Sys.getenv("SGH_PROJ", unset = "C:/Users/rober/OneDrive/Masterkram/R"),
            "scripts", "Finckh_Carmona_2023_PANGAEA_R1.xlsx"),
  format = "file"
)

## ---- foundation: source the scripts, hand down the unit grid -----------------
## Sources 00_setup.R (parameters + null_models_vRN.R), 01_data.R (data load,
## MDL floor, toxic units, unit definitions) and 02_functions.R (all analysis
## functions), then returns the unit grid. `data_files` is listed as a
## dependency purely so the hash of the xlsx gates every refit.
tar_target(
  foundation,
  {
    source(file.path("R", "00_setup.R"))
    source(file.path("R", "01_data.R"))
    source(file.path("R", "02_functions.R"))
    all_units <- c(units, basin_units, section_units, reach_units)
    grid <- expand.grid(unit_i = seq_along(all_units),
                        kind   = c("conc", "tu"),
                        stringsAsFactors = FALSE)
    list(units = all_units, grid = grid, data_file = data_files)
  }
)

## split the grid out so dynamic branching maps over one row per unit x kind
tar_target(unit_grid, foundation$grid)

## ---- unit fits: prep -> glasso+StARS -> modules -> metrics -> cache ----------
## One branch per unit x kind row. Each branch re-sources the three scripts
## (see the bridge caveat in the header): the sourced globals live in the
## branch's evaluation environment, so every branch is self-contained.
## Mirrors the glasso core of run_unit() in 03_mainline.Rmd (the MB/WGCNA
## cross-checks, tables and figures stay in the Rmd pipeline); the returned
## path is the same fitted_*.rds that save_fitted() writes there, tracked via
## format = "file".
tar_target(
  unit_fits,
  {
    source(file.path("R", "00_setup.R"))
    source(file.path("R", "01_data.R"))
    source(file.path("R", "02_functions.R"))

    u  <- foundation$units[[unit_grid$unit_i]]
    pu <- prep_unit(unit_grid$kind, u$site_ids)
    stopifnot(pu$n_comp >= 5, pu$n_sites >= 4)      # same guard as run_unit()

    gl   <- fit_glasso(pu$Xlog)                     # default rep.num = STARS_REPNUM
    memb <- greedy_modules(gl$A)
    m    <- metrics_from_adj(gl$A, memb, pool_n = length(KEEP_GLOBAL[[unit_grid$kind]]))

    ## [repo 15.09.2026] NOTE -- this writes a GLASSO-ONLY fitted artefact
    ## (targets runs the structural core, not the full run_unit() with its MB +
    ## WGCNA branches and tables/figures). S2/S3/S4/S5 read only the glasso
    ## fields, so they are served; the 99_utils figure refresh detects the
    ## missing A_mb/A_wgcna fields and falls back to refitting for such units.
    save_fitted(list(A_glasso = gl$A, mods_glasso = memb,
                     lambda_glasso = gl$lambda, stab_glasso = gl$stab,
                     metrics_glasso = metric_row(m),
                     kind = unit_grid$kind, unit = u,
                     n_sites = pu$n_sites, n_comp = pu$n_comp,
                     rep.num = gl$rep.num),
                unit_grid$kind, u$safe)             # returns the .rds path (invisibly)
  },
  pattern = map(unit_grid),
  format = "file"
)

## ---- run manifest -------------------------------------------------------------
## Depends on unit_fits so it is written last; records the parameter set,
## package versions and git commit for exactly this run (see 02_functions.R).
tar_target(
  manifest,
  {
    ## Each branch is a self-contained environment, so source the foundation
    ## here too (write_run_manifest() lives in 02_functions.R, which computes
    ## KEEP_GLOBAL from 01_data.R's globals at source time -- sourcing
    ## 02_functions.R standalone would fail). Safe: this target depends on
    ## unit_fits, which already requires the data.
    source(file.path("R", "00_setup.R"))
    source(file.path("R", "01_data.R"))
    source(file.path("R", "02_functions.R"))
    m <- write_run_manifest(extra = list(pipeline = "targets",
                                         n_fitted = length(unit_fits)))
    m$timestamp   # write_run_manifest() writes run_manifest.json (or .rds/.txt
                  # without jsonlite) under outroot as a side effect
  }
)

list(data_files, foundation, unit_grid, unit_fits, manifest)
