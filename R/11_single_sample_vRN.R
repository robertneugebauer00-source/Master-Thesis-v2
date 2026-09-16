## =============================================================================
## 11_single_sample_vRN.R -- per-SITE networks: LIONESS-over-glasso and BONOBO
## Robert Neugebauer (vRN). Created 15.09.2026.
##
## WHY THIS FILE EXISTS.
##   The whole pipeline's smallest object is a spatial UNIT (a group of sites) --
##   one glasso network per unit. The SGH test in S5 is therefore a unit-level
##   analogue of Hernandez et al. (2021, ISME J 15:1722), whose ACTUAL design is
##   per-SITE: modularity Q of each site's network vs that site's stress rank,
##   by Spearman. This file builds the per-site networks that let S7 run that
##   test as Hernandez ran it, with n in the hundreds instead of ~21 units.
##
##   WHAT A SINGLE-SITE NETWORK IS -- AND IS NOT (read before interpreting).
##   A PANGAEA site is one grab sample: a single observation vector. You cannot
##   estimate a covariance between two compounds from ONE observation. So a
##   "single-site network" here is NOT "the chemical associations at this site";
##   it is site i's CONTRIBUTION TO / PERTURBATION OF the shared cross-site
##   network. LIONESS makes this explicit (a leave-one-out delta); BONOBO's
##   posterior is the population covariance tilted by site i's rank-1 outer
##   product. Frame it as a per-site readout against the population, not an
##   intrinsic site network -- S7 says exactly this.
##
##   TWO METHODS, TWO CHARACTERS (complementary, not redundant).
##   * LIONESS-over-glasso -> per-site PARTIAL-correlation network. Keeps the
##     conditioning that is the entire reason this project uses glasso, so the
##     pollution/dilution gradient is still conditioned out. This is the method
##     for the per-site SGH test.
##   * BONOBO -> per-site MARGINAL correlation network WITH per-edge uncertainty
##     (a closed-form posterior). Marginal, so it does NOT condition out the
##     gradient -- use it as a cross-check and for per-edge confidence, not as
##     the primary SGH network. bonobo_to_partial() bridges to partials by ridge
##     inversion, with the caveat spelled out at that function.
##
## DEPENDS ON (sourced first, as in every S-notebook):
##   00_setup.R (huge, igraph, STARS_REPNUM), 01_data.R (tu_mat for stress),
##   02_functions.R (prep_unit, fit_glasso, greedy_modules, metrics_from_adj).
## =============================================================================

## null-coalescing helper (base R only gained %||% in 4.4.0; define it if absent)
if (!exists("%||%")) `%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || (length(a) == 1 && is.na(a))) b else a

## ---------------------------------------------------------------------------
## 1. glasso at a FIXED lambda (no StARS) -- the LIONESS base estimator
## ---------------------------------------------------------------------------
## LIONESS is exact only for a linear estimator. glasso is non-linear and
## regularised, so the interpolation is an APPROXIMATION -- and it is only
## defensible if lambda is held FIXED across the full-data fit and every
## leave-one-out fit. Re-selecting lambda per leave-one-out (StARS inside the
## loop) injects selection variance that swamps the single-site signal. So we
## select lambda ONCE on the full unit (fit_glasso -> StARS) and reuse it here.
## @KEY  fixed-lambda glasso partial correlations (the LIONESS base) ----------------------------------
glasso_pcor_at_lambda <- function(Xlog, lambda) {
  h <- huge::huge(as.matrix(Xlog), lambda = lambda, method = "glasso", verbose = FALSE)
  Omega <- as.matrix(h$icov[[1]])            # precision at this single lambda (already L1-sparse)
  d <- sqrt(diag(Omega)); d[d == 0] <- 1
  pcor <- -Omega / outer(d, d); diag(pcor) <- 0
  dimnames(pcor) <- list(colnames(Xlog), colnames(Xlog))
  (pcor + t(pcor)) / 2
}

## ---------------------------------------------------------------------------
## 2. LIONESS over glasso -- one partial-correlation network per site
## ---------------------------------------------------------------------------
## Kuijjer et al. (2019, iScience 14:226). For an estimator f() applied to the
## N-site data and to the N-1 data with site i removed:
##       A_i = N * f(all) - (N-1) * f(all \ i)
## Every A_i shares the aggregate f(all), so the site networks are HIGHLY
## correlated by construction -- S7's SGH regression must use a permutation null
## (permute stress across sites), never parametric p-values. See sgh_persite().
## @KEY  LIONESS single-site partial-correlation networks -----------------------------------------------
lioness_glasso <- function(pu, lambda = NULL, rep.num = STARS_REPNUM,
                           subset_sites = NULL, verbose = TRUE) {
  X <- as.matrix(pu$Xlog)                     # sites x compounds (obs x vars) -- huge wants this
  N <- nrow(X)
  if (N < 4) stop("lioness_glasso: need >= 4 sites")
  if (is.null(lambda)) lambda <- fit_glasso(X, rep.num = rep.num)$lambda   # StARS ONCE, full data
  A_all <- glasso_pcor_at_lambda(X, lambda)
  idx <- if (is.null(subset_sites)) seq_len(N) else which(rownames(X) %in% subset_sites)
  single <- vector("list", length(idx)); names(single) <- rownames(X)[idx]
  for (k in seq_along(idx)) {
    A_loo <- glasso_pcor_at_lambda(X[-idx[k], , drop = FALSE], lambda)
    single[[k]] <- N * A_all - (N - 1) * A_loo
    if (verbose && k %% 10 == 0) message("[lioness] ", k, "/", length(idx), " sites")
  }
  list(agg = A_all, single = single, lambda = lambda, N = N, sites = names(single))
}

## ---------------------------------------------------------------------------
## 3. BONOBO -- faithful port of netZooPy's compute_bonobo()
## ---------------------------------------------------------------------------
## Saha et al. (2024, Genome Res 34:1397). Ported line-for-line from
## netZooPy/netZooPy/bonobo/bonobo.py (verified numerically identical, incl. the
## one trap: numpy's var() is POPULATION variance (ddof=0), so `delta` below uses
## the population variance of the per-compound variances, NOT R's sample var()).
##
## E is variables x observations (COMPOUNDS x SITES) -- i.e. t(pu$Xlog). For
## sample i: prior covariance = cov over all OTHER sites; posterior sample
## covariance = delta * (centred outer product of site i) + (1-delta) * prior;
## coexpression = that covariance scaled to a correlation. delta in (0, 1/3] is
## tuned from the data. Optional sparsification uses the closed-form posterior
## z-score of each covariance entry (inverse-Wishart, posterior d.f. = g + 1/delta).
## @KEY  BONOBO sample-specific coexpression (faithful netZooPy port) ----------------------------------
bonobo_core <- function(E, sample_idx, expression_mean = NULL, delta = NULL,
                        sparsify = FALSE, confidence = 0.05) {
  stopifnot(is.matrix(E)); g <- nrow(E); n <- ncol(E)
  if (is.null(expression_mean)) expression_mean <- rowMeans(E)
  incl <- rep(TRUE, n); incl[sample_idx] <- FALSE
  S <- stats::cov(t(E[, incl, drop = FALSE]))            # g x g prior (leave-one-out), ddof=1 like np.cov
  dvar <- diag(S)
  if (is.null(delta)) {
    var_pop <- mean((dvar - mean(dvar))^2)               # POPULATION var, to match numpy .var()
    delta <- 1 / (3 + 2 * mean(sqrt(dvar)) / var_pop)
  }
  xc <- E[, sample_idx] - expression_mean
  sscov <- delta * outer(xc, xc) + (1 - delta) * S       # sample-specific covariance
  dsd <- sqrt(diag(sscov)); dsd[dsd == 0] <- 1
  coexpr <- sscov / outer(dsd, dsd); diag(coexpr) <- 1   # -> correlation
  pval <- NULL
  if (sparsify) {
    d  <- g + 1 / delta
    a1 <- (d - g + 1) / ((d - g) * (d - g - 3))
    a2 <- (d - g - 1) / ((d - g) * (d - g - 3))
    v  <- diag(sscov)
    z  <- sscov / sqrt(a1 * (sscov * sscov) + a2 * outer(v, v))
    pval <- 2 * (1 - stats::pnorm(abs(z)))
    keep <- abs(z) > stats::qnorm(1 - confidence / 2)
    coexpr <- diag(g) + (1 - diag(g)) * (coexpr * keep)  # sparse off-diagonal, unit diagonal
  }
  dimnames(coexpr) <- list(rownames(E), rownames(E))
  list(coexpr = coexpr, sscov = sscov, delta = delta, pval = pval)
}

## BONOBO for every site of a unit. sparsify=TRUE by default because a DENSE
## correlation matrix has no meaningful modularity -- topology needs the
## z-score-thresholded (sparse) network. Returns marginal correlation networks.
## @KEY  BONOBO across all sites of a unit ------------------------------------------------------------
bonobo_all <- function(pu, sparsify = TRUE, confidence = 0.05, verbose = TRUE) {
  E <- t(as.matrix(pu$Xlog))                 # compounds x sites
  em <- rowMeans(E); N <- ncol(E)
  single <- vector("list", N); names(single) <- colnames(E)
  deltas <- setNames(numeric(N), colnames(E))
  for (i in seq_len(N)) {
    b <- bonobo_core(E, i, expression_mean = em, sparsify = sparsify, confidence = confidence)
    single[[i]] <- b$coexpr; deltas[i] <- b$delta
    if (verbose && i %% 25 == 0) message("[bonobo] ", i, "/", N, " sites")
  }
  list(single = single, delta = deltas, N = N, sites = colnames(E))
}

## OPTIONAL bridge: BONOBO covariance -> partial correlations, by RIDGE inversion.
## CAVEAT -- BONOBO's covariance is not sparse and, when p (compounds) >= n
## (sites), the leave-one-out prior is singular, so a plain inverse does not
## exist. Ridge (add eps to the diagonal) always inverts but the eps is an extra
## knob the method never specified; the partials are therefore only as trustworthy
## as eps. For per-site PARTIAL correlations prefer LIONESS-over-glasso, whose L1
## penalty handles p >= n natively. This is provided for completeness / comparison.
bonobo_to_partial <- function(pu, sample_idx, eps = 1e-2) {
  E <- t(as.matrix(pu$Xlog))
  b <- bonobo_core(E, sample_idx, sparsify = FALSE)
  S <- b$sscov; p <- nrow(S)
  Theta <- solve(S + eps * diag(p))          # ridge-regularised precision
  d <- sqrt(diag(Theta)); d[d == 0] <- 1
  pcor <- -Theta / outer(d, d); diag(pcor) <- 0
  dimnames(pcor) <- dimnames(S)
  (pcor + t(pcor)) / 2
}

## ---------------------------------------------------------------------------
## 4. per-site stress axis -- identical definition to S5
## ---------------------------------------------------------------------------
## S5's stress axis is per-site sum-TU: colSums(tu_mat[present, sites]). Same
## here, so the per-site SGH test in S7 is on the same scale as the unit-level
## one in S5 (no Diazinon special-casing: this repo does not exclude it, and
## consistency with S5 is what matters).
site_stress <- function(site_ids) {
  if (!exists("tu_mat") || !length(tu_mat)) return(setNames(rep(NA_real_, length(site_ids)), site_ids))
  s <- intersect(site_ids, colnames(tu_mat))
  out <- setNames(rep(NA_real_, length(site_ids)), site_ids)
  if (length(s)) out[s] <- colSums(tu_mat[, s, drop = FALSE])
  out
}

## ---------------------------------------------------------------------------
## 5. per-site topology + distance-from-aggregate
## ---------------------------------------------------------------------------
## For each site network A_i: modularity Q (greedy_modules on |A_i|, the repo
## convention), edges, density, mean strength, and how far the site network sits
## from the aggregate (Frobenius norm of the edge difference, and 1 - Pearson
## correlation of the two edge vectors). The distance is a novel per-site
## descriptor: how idiosyncratic is this site's chemical structure vs its basin.
## @KEY  per-site network metrics + idiosyncrasy distance ---------------------------------------------
single_site_metrics <- function(ss, stress = NULL) {
  agg <- ss$agg                                    # NULL for BONOBO (no single aggregate)
  ut  <- if (!is.null(agg)) upper.tri(agg) else NULL
  rows <- lapply(names(ss$single), function(sid) {
    A <- ss$single[[sid]]
    memb <- tryCatch(greedy_modules(A), error = function(e) NULL)
    m <- if (!is.null(memb)) metrics_from_adj(A, memb) else NULL
    dfa_fro <- dfa_cor <- NA_real_
    if (!is.null(agg)) {
      dfa_fro <- sqrt(sum((A[ut] - agg[ut])^2))
      sdA <- stats::sd(A[ut]); sdG <- stats::sd(agg[ut])
      dfa_cor <- if (sdA > 0 && sdG > 0) 1 - stats::cor(A[ut], agg[ut]) else NA_real_
    }
    data.frame(site = sid,
               Q          = if (!is.null(m)) m$modularity else NA_real_,
               edges      = if (!is.null(m)) m$edges      else NA_integer_,
               density    = if (!is.null(m)) m$density    else NA_real_,
               mean_clust = if (!is.null(m)) m$mean_clust else NA_real_,
               strength   = sum(abs(A[upper.tri(A)])),
               dist_fro   = dfa_fro,
               dist_cor   = dfa_cor,
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  if (!is.null(stress)) out$stress <- stress[out$site]
  out
}

## ---------------------------------------------------------------------------
## 6. per-site SGH test -- Spearman Q ~ stress with a PERMUTATION null
## ---------------------------------------------------------------------------
## The per-site networks are not independent (they share the aggregate, LIONESS;
## or the population prior, BONOBO), so parametric Spearman p-values are
## anticonservative. The honest null permutes the stress labels across sites and
## rebuilds the rho distribution. rho < 0 is the SGH-predicted direction
## (modularity falls as stress rises). Also reports a magnitude control: Q is
## mechanically tied to edge count/strength, so a Q-stress slope that vanishes
## once `edges` is partialled out is a density effect, not modular reorganisation.
## @KEY  per-site SGH test with label-permutation null ------------------------------------------------
sgh_persite <- function(metrics_df, B = 4999, seed = 1) {
  d <- metrics_df[is.finite(metrics_df$Q) & is.finite(metrics_df$stress), ]
  if (nrow(d) < 8) return(list(n = nrow(d), note = "too few sites with Q and stress"))
  rho <- suppressWarnings(stats::cor(d$Q, d$stress, method = "spearman"))
  set.seed(seed)
  null <- replicate(B, suppressWarnings(stats::cor(d$Q, sample(d$stress), method = "spearman")))
  p_perm <- (1 + sum(abs(null) >= abs(rho))) / (B + 1)          # two-sided
  ## partial Spearman of Q~stress given edges (rank-based), as a magnitude control
  rr <- data.frame(Q = rank(d$Q), s = rank(d$stress), e = rank(d$edges))
  pr <- tryCatch({
    rqs <- resid(lm(Q ~ e, rr)); rss <- resid(lm(s ~ e, rr)); stats::cor(rqs, rss)
  }, error = function(e) NA_real_)
  list(n = nrow(d), rho = rho, p_perm = p_perm,
       rho_partial_given_edges = pr,
       sgh_direction = if (!is.na(rho) && rho < 0) "consistent (Q falls as stress rises)" else "not consistent",
       null_mean = mean(null), null_sd = stats::sd(null))
}
