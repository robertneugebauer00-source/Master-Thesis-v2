## =============================================================================
## 05_lioness.R -- LIONESS / BONOBO single-sample (per-site) networks
## Robert Neugebauer, MSc thesis -- LIONESS/BONOBO addition 16.09.2026
## Created 16.09.2026 for the S6 supplement (analysis/S6_single_networks.Rmd);
## not split out of the monolithic script -- new code, repo conventions kept.
##
## WHAT THIS FILE IS. The mainline network for a spatial unit is ONE graph
## estimated from ALL of the unit's sites (02_functions.R, fit_glasso()).
## Single-sample network methods instead ask: what does the network look like
## from the perspective of ONE site? Two estimators are implemented:
##
##   LIONESS  (Kuijjer et al. 2019, iScience 14:226): linear extrapolation of
##     the leave-one-out (LOO) estimates, e_q = N*G_all - (N-1)*G_-q, applied
##     here to Pearson correlations (lioness_cor) and to glasso partial
##     correlations at the FIXED aggregate StARS lambda (lioness_glasso).
##   BONOBO   (Saha, Fanfani et al. 2024, Genome Research 35,
##     doi:10.1101/gr.279117.124): Bayesian sample-specific covariance,
##     Sigma_q = delta_q * dx_q dx_q' + (1 - delta_q) * S_-q, replicating the
##     reference implementation netZooPy::compute_bonobo EXACTLY (cross-checked
##     against the Python source, tests/lioness_test.R test 4).
##
## SIDE EFFECTS AT SOURCE TIME: none -- function definitions only. Everything
## that touches repo globals (kind_dir, outroot) is evaluated lazily at call
## time with exists() fallbacks, so this file is sourceable standalone (the
## test suite does exactly that).
## =============================================================================


## =============================================================================
## [S6] ---- closed-form leave-one-out Pearson correlation ----------------------
##
## Brute force would recompute cor(X[-q, ]) for every left-out site q -- an
## O(n^2 p^2) loop of cor() calls. The sufficient statistics of the Pearson
## correlation (column sums, column sums of squares, crossproducts) are
## ADDITIVE over rows, so the LOO statistics are one rank-1 downdate away from
## the full-data ones: O(p^2) per site, O(n p^2) overall.
## =============================================================================

## .loo_cor_one(X, q, S1, S2, Sxy)
##   X   : n sites x p compounds numeric matrix (already log-transformed,
##         MDL/2-floored -- i.e. exactly prep_unit()$Xlog)
##   q   : index of the left-out row (site)
##   S1  : colSums(X), S2 : colSums(X^2), Sxy : crossprod(X)  -- precomputed
##         once by the caller so n LOO rounds cost O(n p^2), not O(n^2 p^2)
##   returns the p x p Pearson correlation matrix of X[-q, ], symmetrised,
##   diagonal 1. Pairs involving a zero-variance column (after the LOO removal)
##   are set to 0 -- cor() would return NA there; 0 keeps the downstream
##   LIONESS extrapolation finite and is the conservative choice (no edge).
.loo_cor_one <- function(X, q, S1, S2, Sxy) {
  n1  <- nrow(X) - 1
  x   <- X[q, ]
  s1  <- S1 - x                    ## column sums without site q
  s2  <- S2 - x^2                  ## column sums of squares without site q
  sxy <- Sxy - tcrossprod(x)       ## crossproducts without site q

  num <- sxy - tcrossprod(s1) / n1             ## centred crossproducts
  den <- s2 - s1^2 / n1                        ## centred sums of squares
  bad <- !is.finite(den) | den <= 0            ## zero-variance guard
  den[bad] <- NA_real_
  r <- num / sqrt(outer(den, den))
  r[is.na(r)] <- 0
  r <- (r + t(r)) / 2                          ## kill fp asymmetry
  diag(r) <- 1
  r
}


## =============================================================================
## [S6] ---- LIONESS on Pearson correlations ------------------------------------
##
## Kuijjer et al. (2019, iScience 14:226): if G_all is a network estimated from
## all N samples and G_-q the same estimate with sample q removed, and the
## estimator is (approximately) LINEAR in the samples, then the contribution of
## sample q is recovered by extrapolation:
##
##     e_q = N * G_all - (N - 1) * G_-q
##
## The mean of the e_q over q equals G_all exactly for statistics LINEAR in
## the samples (e.g. covariance); for Pearson correlation it differs by the
## O(1/n) jackknife bias of the statistic -- both measured and documented in
## tests/lioness_test.R, test 2.
##
## @CAVEAT  LIONESS single-sample networks are NOT correlation matrices: the
## extrapolation amplifies the LOO difference by ~N, so entries leave [-1, 1]
## and the matrices are not PSD. That is a property of the method, not a bug --
## it is precisely what BONOBO (below) was designed to fix. Do not feed
## lioness nets to anything that assumes a valid correlation matrix (e.g.
## eigen-decomposition-based centralities) without thinking; on-skeleton
## summaries (single_site_summary) are safe. `clip = TRUE` truncates entries
## to [-1, 1] for display purposes only -- off by default, because clipping
## destroys the consistency property.
## =============================================================================

## lioness_cor(X, keep, progress, clip)
##   X    : n sites x p compounds matrix (prep_unit()$Xlog)
##   keep : "array" -> nets is a p x p x n array (3rd dim = sites);
##          "list"  -> nets is a named list of n p x p matrices
##   clip : truncate per-site entries to [-1, 1] (display only, see @CAVEAT)
##   returns list(nets, aggregate = cor(X), method = "lioness_pearson")
lioness_cor <- function(X, keep = c("array", "list"), progress = TRUE, clip = FALSE) {
  keep <- match.arg(keep)
  X <- as.matrix(X)
  n <- nrow(X); p <- ncol(X)
  if (n < 4) stop("lioness_cor: need at least 4 sites for a leave-one-out correlation")
  if (is.null(rownames(X))) rownames(X) <- paste0("site_", seq_len(n))
  if (is.null(colnames(X))) colnames(X) <- paste0("V", seq_len(p))

  G_all <- stats::cor(X)
  S1 <- colSums(X); S2 <- colSums(X^2); Sxy <- crossprod(X)

  nets <- array(NA_real_, dim = c(p, p, n),
                dimnames = list(colnames(X), colnames(X), rownames(X)))
  for (q in seq_len(n)) {
    G_loo <- .loo_cor_one(X, q, S1, S2, Sxy)
    e_q   <- n * G_all - (n - 1) * G_loo
    if (clip) e_q <- pmax(pmin(e_q, 1), -1)
    nets[, , q] <- e_q
    if (progress && (q %% max(1, n %/% 10) == 0 || q == n))
      message(sprintf("[S6] lioness_cor: site %d/%d", q, n))
  }

  if (keep == "list")
    nets <- setNames(lapply(seq_len(n), function(q) nets[, , q]), rownames(X))

  list(nets = nets, aggregate = G_all, method = "lioness_pearson")
}


## =============================================================================
## [S6] ---- LIONESS on glasso partial correlations -----------------------------
##
## The mainline network is a glasso PARTIAL-correlation network, so the
## single-sample companion that is directly comparable to it applies the
## LIONESS extrapolation to glasso estimates, not to marginal correlations.
## StARS is NOT re-run per site -- that would be n full stability selections.
## Instead every LOO refit uses the FIXED aggregate lambda (the one StARS
## picked on all sites), so each refit is a single glasso solve: n glasso
## fits at one lambda each, roughly the cost of one mainline fit's lambda
## path. Fixed lambda is also the statistically sensible choice: per-site
## networks then describe the same penalised model, differing only in which
## site is extrapolated.
##
## Edge conversion is EXACTLY fit_glasso's (02_functions.R):
##   pcor = -Omega / outer(sqrt(diag(Omega)), sqrt(diag(Omega))), diag = 0,
##   symmetrised. The glasso solution at a fixed lambda is already sparse, so
##   no extra refit-mask is needed (fit_glasso's mask only re-zeros what StARS
##   deselected along the path).
##
## restrict_to_aggregate = TRUE zeroes every entry whose aggregate edge is
## zero: all per-site networks then live on the SAME StARS-selected skeleton,
## which is what makes their on-skeleton weights comparable across sites (and
## is the analogue of the fixed global node set, one level down).
## =============================================================================

## lioness_glasso(X, lambda, A_agg, restrict_to_aggregate, seed, progress)
##   X      : n sites x p compounds matrix (prep_unit()$Xlog)
##   lambda : FIXED penalty, normally fit_glasso()$lambda / the cached
##            fitted_<kind>_<safe>.rds $lambda_glasso
##   A_agg  : aggregate p x p partial-correlation adjacency; NULL -> fitted
##            here on all sites at the same fixed lambda
##   returns list(nets = p x p x n array, aggregate = A_agg, lambda,
##                restricted, method = "lioness_glasso")
lioness_glasso <- function(X, lambda, A_agg = NULL, restrict_to_aggregate = TRUE,
                           seed = 42, progress = TRUE) {
  if (!requireNamespace("huge", quietly = TRUE))
    stop("lioness_glasso: package 'huge' is required (glasso refits); ",
         "install it or use lioness_cor()/bonobo_cor() instead")
  X <- as.matrix(X)
  n <- nrow(X); p <- ncol(X)
  if (n < 5) stop("lioness_glasso: need at least 5 sites for leave-one-out glasso")
  if (is.null(rownames(X))) rownames(X) <- paste0("site_", seq_len(n))
  if (is.null(colnames(X))) colnames(X) <- paste0("V", seq_len(p))

  ## [vRN-fix pattern 12.08.2026] glasso at a fixed lambda is deterministic,
  ## but keep the seed argument for parity with fit_glasso() and the same
  ## no-RNG-leak guard, so loops over units stay reproducible regardless.
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) stats::runif(1)
  .old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit(assign(".Random.seed", .old_seed, envir = .GlobalEnv), add = TRUE)
  set.seed(seed)

  ## one glasso solve at the fixed lambda -> partial correlations, fit_glasso-style
  .pcor_at_lambda <- function(M) {
    h     <- huge::huge(as.matrix(M), method = "glasso", lambda = lambda,
                        cov.output = TRUE, verbose = FALSE)
    Omega <- as.matrix(h$icov[[1]])          # single lambda -> single precision matrix
    d     <- sqrt(diag(Omega))
    pcor  <- -Omega / outer(d, d)
    diag(pcor) <- 0
    (pcor + t(pcor)) / 2
  }

  if (is.null(A_agg)) A_agg <- .pcor_at_lambda(X)
  dimnames(A_agg) <- list(colnames(X), colnames(X))

  nets <- array(NA_real_, dim = c(p, p, n),
                dimnames = list(colnames(X), colnames(X), rownames(X)))
  for (q in seq_len(n)) {
    A_loo <- .pcor_at_lambda(X[-q, , drop = FALSE])
    e_q   <- n * A_agg - (n - 1) * A_loo
    if (restrict_to_aggregate) e_q[A_agg == 0] <- 0   ## @KEY common skeleton
    diag(e_q) <- 0
    nets[, , q] <- e_q
    if (progress && (q %% max(1, n %/% 10) == 0 || q == n))
      message(sprintf("[S6] lioness_glasso: site %d/%d", q, n))
  }

  list(nets = nets, aggregate = A_agg, lambda = lambda,
       restricted = restrict_to_aggregate, method = "lioness_glasso")
}


## =============================================================================
## [S6] ---- BONOBO: Bayesian single-sample correlation networks -----------------
##
## Saha, Fanfani et al. (2024, Genome Research, doi:10.1101/gr.279117.124).
## BONOBO models each sample's covariance as a Bayesian posterior mean between
## the leave-one-out covariance S_-q (the prior, everything the other sites
## say) and the sample's own outer-product deviation dx_q dx_q' (the
## likelihood, what only this site says):
##
##     Sigma_q = delta_q * dx_q dx_q' + (1 - delta_q) * S_-q
##
## with the posterior weight tuned from the data,
##
##     delta_q = 1 / (3 + 2 * mean(sqrt(diag(S_-q))) / var(diag(S_-q)))
##
## [S6] PROVENANCE NOTE: this is the netZooPy-CALIBRATED form of delta
## (netZooPy::compute_bonobo, bonobo.py, including numpy's population variance,
## ddof = 0, in the denominator). The paper's eq. 8 differs slightly; the
## reference implementation is replicated here as the default so results are
## comparable to published netZooPy output -- tests/lioness_test.R test 4
## cross-checks against the Python original to 1e-6.
##
## Unlike LIONESS nets, BONOBO networks are proper correlation matrices:
## Sigma_q is a convex combination of two PSD matrices, hence PSD, so after
## scaling, every per-site network is PSD with entries in [-1, 1] and diagonal
## 1. That is the property that makes BONOBO the safer choice for any
## downstream score that assumes a valid correlation matrix.
##
## Sparsification (optional, also exactly netZooPy): the approximate sampling
## variance of the Sigma_q entries under the posterior is
##     sd_jk = sqrt(a1 * Sigma_q^2 + a2 * outer(diag(Sigma_q), diag(Sigma_q)))
##     d = g + 1/delta_q,  a1 = (d-g+1)/((d-g)(d-g-3)),  a2 = (d-g-1)/((d-g)(d-g-3))
## giving z = Sigma_q / sd_jk, two-sided normal p-values, and a sparse network
## keeping only entries with |z| > qnorm(1 - alpha/2) (diagonal kept at 1).
## =============================================================================

## bonobo_cor(X, delta, alpha, progress, save_pvals, compute_sparse)
##   X             : n sites x p compounds matrix (prep_unit()$Xlog)
##   delta         : NULL -> tuned per site (netZooPy form, see above); or a
##                   single numeric in (0, 1) used for every site
##   alpha         : two-sided confidence level for the sparse mask
##   save_pvals    : also return the p x p x n array of z-test p-values
##                   (memory: one more full array -- off by default)
##   compute_sparse: also return the masked sparse per-site networks
##   returns list(nets = p x p x n dense correlations, deltas = named numeric,
##                pvals = array or NULL, sparse = array or NULL, alpha,
##                method = "bonobo")
bonobo_cor <- function(X, delta = NULL, alpha = 0.05, progress = TRUE,
                       save_pvals = FALSE, compute_sparse = TRUE) {
  X <- as.matrix(X)
  n <- nrow(X); p <- ncol(X)
  if (n < 4) stop("bonobo_cor: need at least 4 sites for a leave-one-out covariance")
  if (is.null(rownames(X))) rownames(X) <- paste0("site_", seq_len(n))
  if (is.null(colnames(X))) colnames(X) <- paste0("V", seq_len(p))
  if (!is.null(delta) && (length(delta) != 1 || !is.finite(delta) || delta <= 0 || delta >= 1))
    stop("bonobo_cor: delta must be NULL or a single numeric in (0, 1)")

  ## netZooPy works genes x samples; mirror that orientation internally
  E    <- t(X)                          # p compounds x n sites
  xbar <- rowMeans(E)                   # full-data mean (over ALL sites), as netZooPy

  dn   <- list(colnames(X), colnames(X), rownames(X))
  nets   <- array(NA_real_, dim = c(p, p, n), dimnames = dn)
  pvals  <- if (save_pvals)     array(NA_real_, dim = c(p, p, n), dimnames = dn) else NULL
  sparse <- if (compute_sparse) array(NA_real_, dim = c(p, p, n), dimnames = dn) else NULL
  deltas <- setNames(rep(NA_real_, n), rownames(X))
  thr    <- stats::qnorm(1 - alpha / 2)

  for (q in seq_len(n)) {
    ## [S6] orientation trap: np.cov treats ROWS as variables, R's cov() treats
    ## COLUMNS as variables -- so this is cov() on the sites-x-compounds side,
    ## not cov(E[, -q]). ddof = 1 in both, as np.cov's default.
    S_q <- stats::cov(X[-q, , drop = FALSE])
    dg  <- diag(S_q)

    if (is.null(delta)) {
      ## [S6] netZooPy form: numpy's var() is the POPULATION variance (ddof=0)
      v_pop <- mean((dg - mean(dg))^2)
      delta_q <- if (is.finite(v_pop) && v_pop > 0)
        1 / (3 + 2 * mean(sqrt(dg)) / v_pop) else 1 / 3   ## flat fallback, documented
    } else delta_q <- delta
    deltas[q] <- delta_q

    dx      <- E[, q] - xbar
    Sigma_q <- delta_q * tcrossprod(dx) + (1 - delta_q) * S_q

    ## covariance -> correlation; zero diagonal replaced by 1 (netZooPy guard)
    sd_q <- sqrt(diag(Sigma_q))
    sd_q[sd_q == 0 | !is.finite(sd_q)] <- 1
    R_q  <- Sigma_q / outer(sd_q, sd_q)
    R_q  <- (R_q + t(R_q)) / 2
    nets[, , q] <- R_q

    if (compute_sparse || save_pvals) {
      g  <- p
      d  <- g + 1 / delta_q
      a1 <- (d - g + 1) / ((d - g) * (d - g - 3))
      a2 <- (d - g - 1) / ((d - g) * (d - g - 3))
      sd_jk <- sqrt(a1 * Sigma_q^2 + a2 * tcrossprod(diag(Sigma_q)))
      z     <- Sigma_q / sd_jk
      pv    <- 2 * (1 - stats::pnorm(abs(z)))
      if (save_pvals) pvals[, , q] <- pv
      if (compute_sparse) {
        mask <- abs(z) > thr
        mask[is.na(mask)] <- FALSE        ## numpy: nan > x is False
        sp  <- R_q * mask
        diag(sp) <- 1                     ## netZooPy: eye + (1-eye) * masked dense
        sparse[, , q] <- sp
      }
    }
    if (progress && (q %% max(1, n %/% 10) == 0 || q == n))
      message(sprintf("[S6] bonobo_cor: site %d/%d", q, n))
  }

  list(nets = nets, deltas = deltas, pvals = pvals, sparse = sparse,
       alpha = alpha, method = "bonobo")
}


## =============================================================================
## [S6] ---- per-site summary scores --------------------------------------------
##
## igraph-free on purpose: the per-site objects are what feeds downstream
## regressions (site-level network score ~ site-level stress, the SGH readout
## in S6), and that needs a tidy data.frame, not a graph object.
## =============================================================================

## single_site_summary(nets, A_agg, membership)
##   nets       : p x p x n array (or named list of n p x p matrices), as
##                returned in $nets by lioness_cor / lioness_glasso / bonobo_cor
##   A_agg      : optional aggregate adjacency; the SKELETON is A_agg != 0
##   membership : optional named module vector (e.g. fit$mods_glasso)
##   returns a data.frame, one row per site (row names = site ids):
##     mean_weight / sd_weight      over ALL off-diagonal entries
##     n_edges_on_skeleton          skeleton edges with non-zero weight at this
##                                  site (= skeleton size for dense nets)
##     mean_weight_on_skeleton      mean per-site weight over skeleton edges
##     strength_skeleton            sum of |weight| over skeleton edges
##     within_module_mean           mean weight over skeleton edges whose both
##                                  endpoints share a module (needs membership)
single_site_summary <- function(nets, A_agg = NULL, membership = NULL) {
  if (is.list(nets) && !is.array(nets)) {
    nn <- names(nets); if (is.null(nn)) nn <- paste0("site_", seq_along(nets))
    p  <- ncol(nets[[1]])
    nets <- array(unlist(nets, use.names = FALSE), dim = c(p, p, length(nets)),
                  dimnames = list(colnames(nets[[1]]), colnames(nets[[1]]), nn))
  }
  p <- dim(nets)[1]; n <- dim(nets)[3]
  sites <- dimnames(nets)[[3]]; if (is.null(sites)) sites <- paste0("site_", seq_len(n))
  comps <- dimnames(nets)[[1]]; if (is.null(comps)) comps <- paste0("V", seq_len(p))

  ut   <- upper.tri(nets[, , 1])
  skel <- if (!is.null(A_agg)) (A_agg != 0) & ut else ut
  n_skel <- sum(skel)

  within <- NULL
  if (!is.null(membership)) {
    mm <- membership[comps]
    if (any(is.na(mm))) {
      ## compounds without a module assignment each form their own singleton,
      ## so they can never contribute a within-module edge
      mm[is.na(mm)] <- max(mm, na.rm = TRUE) + seq_len(sum(is.na(mm)))
    }
    within <- outer(mm, mm, `==`) & skel
  }

  rows <- lapply(seq_len(n), function(q) {
    W  <- nets[, , q]
    w  <- W[ut]
    ws <- W[skel]
    data.frame(
      site                   = sites[q],
      mean_weight            = mean(w),
      sd_weight              = stats::sd(w),
      n_edges_on_skeleton    = if (!is.null(A_agg)) sum(ws != 0) else NA_integer_,
      mean_weight_on_skeleton = if (!is.null(A_agg) && n_skel > 0) mean(ws) else NA_real_,
      strength_skeleton      = if (!is.null(A_agg) && n_skel > 0) sum(abs(ws)) else NA_real_,
      within_module_mean     = if (!is.null(within) && sum(within) > 0) mean(W[within]) else NA_real_,
      stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- out$site
  out
}


## =============================================================================
## [S6] ---- caching helpers (fitted_<kind>_<safe>.rds conventions) -------------
##
## Mirrors fitted_path()/load_fitted()/save_fitted() in 02_functions.R, one
## .rds per unit per data kind per single-sample METHOD. kind_dir() and
## outroot are looked up LAZILY (call time) with exists() fallbacks, so this
## file stays sourceable without 00-02 loaded (tests do exactly that).
## =============================================================================

## single_net_path(kind, safe, method, root) -> path to single_<method>_<kind>_<safe>.rds
single_net_path <- function(kind, safe, method,
                            root = if (exists("outroot")) get("outroot") else "outputs_vRN") {
  kd <- if (exists("kind_dir", mode = "function")) kind_dir(kind) else
        if (kind == "conc") "concentration" else "toxic_unit"
  file.path(root, kd, safe, sprintf("single_%s_%s_%s.rds", method, kind, safe))
}

## save_single_nets(obj, kind, safe, method, root) -- obj is the list returned
## by lioness_cor()/lioness_glasso()/bonobo_cor(); returns the path invisibly.
save_single_nets <- function(obj, kind, safe, method,
                             root = if (exists("outroot")) get("outroot") else "outputs_vRN") {
  f <- single_net_path(kind, safe, method, root)
  dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
  saveRDS(obj, f)
  message(sprintf("[S6] saved %s single-sample nets: %s", method, f))
  invisible(f)
}

## load_single_nets(kind, safe, method, root) -- cached list, or NULL if absent.
load_single_nets <- function(kind, safe, method,
                             root = if (exists("outroot")) get("outroot") else "outputs_vRN") {
  f <- single_net_path(kind, safe, method, root)
  if (file.exists(f)) readRDS(f) else NULL
}
