## =============================================================================
## tests/smoke_test.R -- end-to-end smoke test on SYNTHETIC data
## [repo 15.09.2026]
##
## PURPOSE. Prove that the core pipeline -- prep_unit() -> fit_glasso() (glasso
## + StARS) -> greedy_modules() -> null_test_multi() -- recovers planted
## modular structure, and that fit_glasso() keeps its promise not to leak its
## internal seed into the global RNG. Runs WITHOUT the private PANGAEA data and
## WITHOUT the heavy dependencies: only igraph and huge are required. WGCNA,
## readxl, magick, uwot and openxlsx are deliberately NOT needed (the WGCNA
## branch fit_wgcna() and all plotting functions are skipped -- nothing here
## calls them), and 00_setup.R is deliberately NOT sourced (it loads WGCNA et
## al.); instead the handful of globals that 02_functions.R touches at SOURCE
## TIME are stubbed below on synthetic data.
##
## RUN. From the repo root:
##   Rscript tests/smoke_test.R
## (Any working directory works; the script walks up from its own location
## until it finds run_all.R.) Exit status 0 and "SMOKE TEST PASSED" on success,
## exit status 1 otherwise.
##
## WHAT 02_functions.R EVALUATES AT SOURCE TIME (checked against the file):
##   - REF_SITES <- sites_meta$site_id[sites_meta$river_basin %in% big_basins]
##   - KEEP_GLOBAL <- list(conc = global_nodeset("conc"), tu = global_nodeset("tu"))
##     which needs conc_mat, tu_mat and PREV_THRESH, plus the message() after it.
##   - RISK_RING (a plain constant).
## Everything else (ug_palette/ug_lookup in the plotting functions, tu_mat and
## RQ_THRESH in risk_driver_table(), outroot in save_fitted()/load_fitted(),
## STARS_REPNUM as a fit_glasso() default) is inside function bodies or default
## arguments and is only evaluated at CALL time -- we stub what we call
## (PREV_THRESH, NODESET_SCOPE, MIN_DETECT_UNIT, mdl_half_ngL, mdl_half_tu) and
## never call the rest (we pass rep.num explicitly, so STARS_REPNUM is not
## needed either).
## =============================================================================

## ---- locate the repo root ---------------------------------------------------
## Walk up from the directory containing this script (from Rscript's --file
## argument; fall back to getwd() for interactive source()) until run_all.R.
.args     <- commandArgs(trailingOnly = FALSE)
.file_arg <- sub("^--file=", "", grep("^--file=", .args, value = TRUE))
.start    <- if (length(.file_arg)) dirname(normalizePath(.file_arg[1], winslash = "/"))
             else getwd()
ROOT <- .start
repeat {
  if (file.exists(file.path(ROOT, "run_all.R"))) break
  .parent <- dirname(ROOT)
  if (identical(.parent, ROOT)) stop("smoke_test.R: could not locate repo root (run_all.R not found above ", .start, ")")
  ROOT <- .parent
}
cat("repo root:", ROOT, "\n")

## ---- tiny check harness ------------------------------------------------------
.failures <- 0L
check <- function(ok, msg) {
  ok <- isTRUE(ok)
  cat(sprintf("[%s] %s\n", if (ok) "PASS" else "FAIL", msg))
  if (!ok) .failures <<- .failures + 1L
  invisible(ok)
}

## ---- synthetic data with planted modules ------------------------------------
## 30 sites x 40 compounds. Compounds 1-30 form three planted modules of 10:
## within a module, compounds share a latent log-scale factor. All compounds
## additionally share one WEAKER global pollution-gradient factor. Compounds
## 31-40 are independent apart from that global gradient. Exponentiate to the
## concentration scale, then threshold the lowest ~20% of all values to 0 to
## mimic true non-detects (left censoring), exactly what prep_unit()'s MDL/2
## floor is built for.
set.seed(1)
n_sites <- 30L
n_comp  <- 40L
site_ids <- sprintf("EUS_%03d", seq_len(n_sites))
comp_ids <- sprintf("Compound_%02d", seq_len(n_comp))

g_global <- rnorm(n_sites)                      # global pollution gradient (weak loadings)
Fmod     <- cbind(rnorm(n_sites), rnorm(n_sites), rnorm(n_sites))   # one factor per planted module
truth    <- setNames(rep(1:3, each = 10), comp_ids[1:30])           # planted module membership

X <- matrix(NA_real_, n_sites, n_comp, dimnames = list(site_ids, comp_ids))
for (j in 1:30)
  X[, j] <- 2 + 0.9 * Fmod[, truth[j]] + 0.3 * g_global + rnorm(n_sites, sd = 0.35)
for (j in 31:40)
  X[, j] <- 2 + 0.3 * g_global + rnorm(n_sites, sd = 0.50)

conc <- 10^X
conc[conc <= quantile(conc, 0.20)] <- 0         # ~20% true non-detects
conc_mat <- t(conc)                             # compounds x sites, as in 01_data.R

## ---- stub the globals 02_functions.R needs (see header) ----------------------
## [repo 15.09.2026] minimal stand-ins for the 01_data.R globals.
sites_meta <- data.frame(site_id     = site_ids,
                         river_basin = "Synthetic",
                         lat         = 50 + seq_len(n_sites) / 100,
                         lon         = 8  + seq_len(n_sites) / 100,
                         country     = "XX",
                         stringsAsFactors = FALSE)
big_basins <- "Synthetic"

## TU matrix: same detections, rescaled per compound (TU = MEC / EC10 in the
## real data; here a runif factor per compound plays the EC10 role).
tu_scale <- setNames(runif(n_comp, 1e-6, 1e-3), comp_ids)
tu_mat   <- sweep(conc_mat, 1, tu_scale, "*")

mdl_half_ngL <- setNames(rep(0.5, n_comp), comp_ids)   # one MDL/2 floor per compound
mdl_half_tu  <- mdl_half_ngL * tu_scale                # same floor on the TU scale

PREV_THRESH     <- 0.10
NODESET_SCOPE   <- "global"
MIN_DETECT_UNIT <- 2L

suppressPackageStartupMessages({
  library(igraph)
  library(huge)
})
source(file.path(ROOT, "R", "02_functions.R"))
source(file.path(ROOT, "R", "null_models_vRN.R"))

## ---- 1. prep_unit ------------------------------------------------------------
pu <- prep_unit("conc", site_ids)
check(pu$n_comp >= 30, sprintf("prep_unit: n_comp >= 30 (got %d)", pu$n_comp))
check(pu$n_sites == n_sites, sprintf("prep_unit: n_sites == %d (got %d)", n_sites, pu$n_sites))

## ---- 2 + 3. fit_glasso: valid adjacency, no RNG leak --------------------------
## Seed hygiene: fit_glasso() promises to restore .Random.seed on exit
## (vRN-fix 12.08.2026). Record the global RNG state around the call and
## assert it is bit-identical afterwards.
set.seed(20260915)                              # make sure .Random.seed exists
seed_before <- .Random.seed
gl <- fit_glasso(pu$Xlog, rep.num = 20)         # rep.num = 20 for speed (StARS depth only)
seed_after <- .Random.seed
A <- gl$A

check(is.matrix(A) && nrow(A) == ncol(A) && isTRUE(all.equal(A, t(A), tolerance = 1e-10)),
      "fit_glasso: A is square and symmetric")
check(all(diag(A) == 0), "fit_glasso: diagonal of A is zero")
check(identical(seed_before, seed_after), "fit_glasso: .Random.seed unchanged (no RNG leak)")

## ---- 4. module recovery -------------------------------------------------------
memb <- greedy_modules(A)
deg  <- setNames(colSums(A != 0), colnames(A))
## Compare on the planted compounds that actually have degree > 0 in the fitted
## graph: isolates carry no co-occurrence information, so a planted node that
## StARS left unconnected cannot be assigned correctly and is excluded here.
## (With these signal strengths essentially all 30 planted nodes have edges;
## the < 15 fallback keeps the test meaningful even on an unlucky StARS draw.)
ev <- names(truth)[deg[names(truth)] > 0]
if (length(ev) < 15) ev <- names(truth)
obs_ari <- ari(memb[ev], truth[ev])
check(is.finite(obs_ari) && obs_ari >= 0.6,
      sprintf("greedy_modules: ARI vs planted truth >= 0.6 on %d connected planted compounds (got %.3f)",
              length(ev), obs_ari))

## ---- 5. null models ------------------------------------------------------------
nt <- null_test_multi(A, memb, N = 50)          # 3 nulls x 2 metrics
check(is.data.frame(nt) && nrow(nt) == 6,
      sprintf("null_test_multi: 6 rows (3 nulls x 2 metrics) (got %d)", nrow(nt)))
cfg <- nt[nt$nullModel == "config" & nt$metric == "Modularity", ]
check(nrow(cfg) == 1 && is.finite(cfg$observed) && cfg$observed > cfg$null_mean,
      sprintf("null_test_multi: observed Q (%.3f) exceeds config null mean (%.3f)",
              cfg$observed, cfg$null_mean))

## ---- 6. signed-edge check -------------------------------------------------------
sgn <- signed_edge_check(A)
check(is.data.frame(sgn) && nrow(sgn) == 1 && "pct_negative" %in% names(sgn),
      "signed_edge_check: one-row data.frame with column pct_negative")

## ---- 7. metric row ---------------------------------------------------------------
mr <- metric_row(metrics_from_adj(A, memb))
check(is.numeric(mr) && !is.null(names(mr)) &&
        all(c("Modularity Q", "Nodes (pool)") %in% names(mr)),
      "metric_row: named numeric containing 'Modularity Q' and 'Nodes (pool)'")

## ---- verdict ----------------------------------------------------------------------
if (.failures > 0L) {
  cat(.failures, "check(s) FAILED\n")
  quit(status = 1)
}
cat("SMOKE TEST PASSED\n")
