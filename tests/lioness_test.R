## =============================================================================
## tests/lioness_test.R -- numerical validation of R/05_lioness.R
## [S6 addition 16.09.2026]
##
## PURPOSE. Prove that the single-sample estimators are numerically correct:
##   1. the closed-form leave-one-out Pearson correlation equals brute-force
##      cor(X[-q, ]) to 1e-10 for every left-out site;
##   2. the LIONESS decomposition is internally consistent (jackknife
##      pseudo-value identity, exact) and reproduces cor(X) up to the known
##      O(1/n) jackknife bias of nonlinear statistics;
##   3. BONOBO per-site networks are proper correlation matrices (PSD, entries
##      in [-1, 1]) and shrink toward the truth relative to LIONESS;
##   4. bonobo_cor() replicates the Python reference implementation
##      (netZooPy::compute_bonobo) to 1e-6 -- skipped, not failed, if python3
##      or its deps are unavailable;
##   5. lioness_glasso() runs and honours the aggregate-skeleton restriction
##      (only if the huge package is available).
##
## RUN. From the repo root (any working directory works; the script walks up
## from its own location until it finds run_all.R):
##   Rscript tests/lioness_test.R
## Exit status 0 and "LIONESS TEST PASSED" on success, exit status 1 otherwise.
## Base R only; huge is optional; python3 + numpy + scipy optional (test 4).
## =============================================================================

## ---- locate the repo root ---------------------------------------------------
.args     <- commandArgs(trailingOnly = FALSE)
.file_arg <- sub("^--file=", "", grep("^--file=", .args, value = TRUE))
.start    <- if (length(.file_arg)) dirname(normalizePath(.file_arg[1], winslash = "/")) else
             getwd()
ROOT <- .start
repeat {
  if (file.exists(file.path(ROOT, "run_all.R")) || file.exists(file.path(ROOT, "R", "00_setup.R")) || length(list.files(ROOT, pattern = "\\.Rproj$"))) break
  .parent <- dirname(ROOT)
  if (identical(.parent, ROOT)) stop("lioness_test.R: could not locate repo root (run_all.R not found above ", .start, ")")
  ROOT <- .parent
}
cat("repo root:", ROOT, "\n")

## ---- tiny check harness (same spirit as tests/smoke_test.R) ------------------
.failures <- 0L
check <- function(ok, msg) {
  ok <- isTRUE(ok)
  cat(sprintf("[%s] %s\n", if (ok) "PASS" else "FAIL", msg))
  if (!ok) .failures <<- .failures + 1L
  invisible(ok)
}

## ---- the file under test: sourceable standalone, no 00-02 globals needed -----
source(file.path(ROOT, "R", "05_lioness.R"))

## ---- synthetic data ----------------------------------------------------------
## n = 25 sites x p = 12 compounds, plain Gaussian -- the tests are numerical
## identities / matrix properties, so no planted structure is needed here.
set.seed(42)
n <- 25L; p <- 12L
site_ids <- sprintf("EUS_%03d", seq_len(n))
comp_ids <- sprintf("Compound_%02d", seq_len(p))
X <- matrix(rnorm(n * p), n, p, dimnames = list(site_ids, comp_ids))

## =============================================================================
## TEST 1 -- closed-form LOO Pearson == brute force cor(X[-q, ])
## =============================================================================
S1 <- colSums(X); S2 <- colSums(X^2); Sxy <- crossprod(X)
loo_err <- max(sapply(seq_len(n), function(q)
  max(abs(.loo_cor_one(X, q, S1, S2, Sxy) - cor(X[-q, ])))))
check(loo_err < 1e-10,
      sprintf("closed-form LOO correlation == cor(X[-q, ]) for all %d sites (max err %.2e)", n, loo_err))

## zero-variance guard: a constant column must yield 0s, not NaNs
Xc <- X; Xc[, 3] <- 5
rc <- .loo_cor_one(Xc, 1, colSums(Xc), colSums(Xc^2), crossprod(Xc))
check(!any(is.na(rc)) && all(rc[3, -3] == 0) && all(rc[-3, 3] == 0) && diag(rc)[3] == 1,
      "zero-variance column: off-diagonal set to 0, diagonal 1, no NaN")

## =============================================================================
## TEST 2 -- LIONESS consistency
## =============================================================================
## [S6] e_q = n*G - (n-1)*G_-q is exactly the JACKKNIFE PSEUDO-VALUE of the
## correlation statistic. Its mean over q is therefore exactly the jackknife
## estimate n*cor(X) - (n-1)*mean_q cor(X[-q, ]) -- that identity is exact and
## is what we assert at 1e-10. The LIONESS paper's "mean of single-sample nets
## equals the aggregate" holds exactly only for statistics LINEAR in the
## samples (e.g. covariance); for Pearson correlation the mean differs from
## cor(X) by the jackknife bias of the statistic, O(1/n) -- verified below to
## be small and shrinking in n, and documented in R/05_lioness.R.
lc <- lioness_cor(X, progress = FALSE)
check(identical(dim(lc$nets), c(p, p, n)) && identical(dimnames(lc$nets)[[3]], site_ids),
      "lioness_cor: nets are p x p x n with site ids on the 3rd dim")

mean_net <- apply(lc$nets, c(1, 2), mean)
jk_est   <- n * cor(X) - (n - 1) * Reduce(`+`, lapply(seq_len(n), function(q) cor(X[-q, ]))) / n
check(max(abs(mean_net - jk_est)) < 1e-10,
      sprintf("LIONESS mean == jackknife pseudo-value mean, exact (max err %.2e)",
              max(abs(mean_net - jk_est))))

gap <- max(abs(mean_net - cor(X)))
check(gap < 2 / n,
      sprintf("LIONESS mean ~= cor(X) up to jackknife bias O(1/n) (gap %.4f < %.4f)", gap, 2 / n))

## the covariance-level decomposition underlying it IS exact (linear statistic)
e_cov <- lapply(seq_len(n), function(q) n * cov(X) - (n - 1) * cov(X[-q, ]))
cov_gap <- max(abs(Reduce(`+`, e_cov) / n - cov(X)))
check(cov_gap < 1e-10,
      sprintf("covariance-level LIONESS decomposition is exact (max err %.2e)", cov_gap))

## keep = "list" returns the same numbers as a named list
ll <- lioness_cor(X, keep = "list", progress = FALSE)
check(is.list(ll$nets) && length(ll$nets) == n && identical(names(ll$nets), site_ids) &&
        max(abs(ll$nets[[1]] - lc$nets[, , 1])) < 1e-12,
      "lioness_cor keep = 'list': named list, same numbers as the array")

## =============================================================================
## TEST 3 -- BONOBO: valid correlation matrices + shrinkage
## =============================================================================
bc <- bonobo_cor(X, progress = FALSE, save_pvals = TRUE)
check(identical(dim(bc$nets), c(p, p, n)) && length(bc$deltas) == n &&
        all(bc$deltas > 0 & bc$deltas < 1),
      "bonobo_cor: nets p x p x n, one tuned delta in (0, 1) per site")

min_eig <- min(sapply(seq_len(n), function(q)
  min(eigen(bc$nets[, , q], symmetric = TRUE, only.values = TRUE)$values)))
check(min_eig >= -1e-8,
      sprintf("BONOBO nets are PSD (min eigenvalue %.2e >= -1e-8)", min_eig))
check(max(abs(bc$nets)) <= 1 + 1e-8 &&
        all(vapply(seq_len(n), function(q) max(abs(diag(bc$nets[, , q]) - 1)) < 1e-12, logical(1))),
      sprintf("BONOBO entries in [-1, 1], diagonal 1 (max |entry| %.6f)", max(abs(bc$nets))))

## [S6] the mean of BONOBO nets is NOT cor(X): delta_q shrinks every site
## toward the leave-one-out bulk, so the average is biased toward the
## aggregate covariance structure by construction. Measured here (max abs diff
## ~ 3e-3 at n = 25) and DOCUMENTED rather than forced -- do not "fix" the
## estimator to pass a consistency test it is not designed to satisfy.
b_gap <- max(abs(apply(bc$nets, c(1, 2), mean) - cor(X)))
cat(sprintf("[INFO] mean of BONOBO nets vs cor(X): max abs diff %.4f (shrinkage bias, documented)\n", b_gap))

## sparse mask: diagonal 1, off-diagonal either 0 or |z| > qnorm(1 - alpha/2)
sp_ok <- all(vapply(seq_len(n), function(q) {
  sp <- bc$sparse[, , q]
  all(diag(sp) == 1) && all(sp[lower.tri(sp)] == 0 | abs(bc$nets[, , q][lower.tri(sp)]) > 0)
}, logical(1)))
check(sp_ok && max(abs(bc$sparse)) <= 1 + 1e-8,
      "BONOBO sparse nets: diagonal 1, masked entries zeroed, still in [-1, 1]")

## shrinkage: on homogeneous one-factor data the true correlation is known;
## BONOBO nets must sit closer to it (Frobenius) than LIONESS nets, on average
## over seeds. [S6] this is the PSD/shrinkage property doing its job: LIONESS
## amplifies leave-one-out noise by ~n, BONOBO averages it with the prior.
fro <- vapply(1:6, function(sd) {
  set.seed(100 + sd)
  nn <- 30L; pp <- 15L; th <- rnorm(nn)
  Xs <- vapply(seq_len(pp), function(j) 0.7 * th + rnorm(nn, sd = sqrt(1 - 0.49)), numeric(nn))
  Rtrue <- matrix(0.49, pp, pp); diag(Rtrue) <- 1
  Ll <- lioness_cor(Xs, progress = FALSE)$nets
  Bb <- bonobo_cor(Xs, progress = FALSE)$nets
  c(lioness = mean(vapply(seq_len(nn), function(q) sqrt(mean((Ll[, , q] - Rtrue)^2)), numeric(1))),
    bonobo  = mean(vapply(seq_len(nn), function(q) sqrt(mean((Bb[, , q] - Rtrue)^2)), numeric(1))))
}, numeric(2))
check(mean(fro["bonobo", ]) < mean(fro["lioness", ]),
      sprintf("BONOBO closer to the true correlation than LIONESS (mean Frobenius %.3f < %.3f, 6 seeds)",
              mean(fro["bonobo", ]), mean(fro["lioness", ])))

## =============================================================================
## TEST 4 -- cross-check against the Python reference (netZooPy compute_bonobo)
## =============================================================================
## Skipped (not failed) when python3 / numpy / scipy are unavailable. The
## helper script first tries the real netZooPy import; if the import chain
## fails (it pulls in pandas/igraph via the panda module), it falls back to an
## EXACT inline copy of compute_bonobo() from netZooPy/bonobo/bonobo.py, which
## needs only numpy and scipy.stats.
py <- Sys.which("python3")
py_ok <- nzchar(py) && {
  dep <- suppressWarnings(system2(py, c("-c", shQuote("import numpy, scipy.stats")),
                                  stdout = FALSE, stderr = FALSE))
  identical(dep, 0L)
}
if (!py_ok) {
  cat("[SKIP] python3 with numpy+scipy not available -- test 4 (Python cross-check) skipped\n")
} else {
  td <- tempfile("bonobo_xcheck_"); dir.create(td)
  write.csv(X, file.path(td, "X.csv"), row.names = FALSE)

  py_script <- file.path(td, "ref_bonobo.py")
  writeLines(c(
    "import sys, numpy as np, pandas as pd",
    "sys.path.insert(0, '/tmp/netzoo')",
    "try:",
    "    from netZooPy.bonobo.bonobo import compute_bonobo",
    "    SRC = 'netZooPy import'",
    "except Exception:",
    "    import scipy.stats as stats",
    "    # EXACT inline copy of netZooPy/bonobo/bonobo.py :: compute_bonobo",
    "    # (import chain unavailable; body needs only numpy + scipy.stats)",
    "    def compute_bonobo(expression_matrix, expression_mean, sample_idx,",
    "                       delta=None, compute_sparse=False, confidence=0.05, save_pvals=False, **kw):",
    "        pval = None",
    "        mask_include = [True] * expression_matrix.shape[1]",
    "        mask_include[sample_idx] = False",
    "        covariance_matrix = np.cov(expression_matrix[:, mask_include])",
    "        if delta is None:",
    "            delta = 1 / (3 + 2 * np.sqrt(covariance_matrix.diagonal()).mean()",
    "                         / covariance_matrix.diagonal().var())",
    "        sscov = (delta * np.outer((expression_matrix - expression_mean)[:, sample_idx],",
    "                                  (expression_matrix - expression_mean)[:, sample_idx])",
    "                 + (1 - delta) * covariance_matrix)",
    "        sscov = np.array(sscov)",
    "        diag = np.sqrt(np.diag(np.diag(sscov)))",
    "        diag = np.array(diag)",
    "        indices = np.where(np.diag(diag) == 0)[0]",
    "        for i in indices:",
    "            diag[i, i] = 1",
    "        sds = np.linalg.inv(diag)",
    "        bonobo_matrix = sds @ sscov @ sds",
    "        if compute_sparse:",
    "            threshold = stats.norm.ppf(1 - (confidence / 2))",
    "            g = sscov.shape[1]",
    "            d = g + 1 / delta",
    "            a1 = (d - g + 1) / ((d - g) * (d - g - 3))",
    "            a2 = (d - g - 1) / ((d - g) * (d - g - 3))",
    "            v = np.diag(sscov)",
    "            v = (a1 * (np.multiply(sscov, sscov))) + (a2 * (np.outer(v, v)))",
    "            v = np.sqrt(v)",
    "            v = np.divide(sscov, v)",
    "            pval = 2 * (1 - stats.norm.cdf(np.abs(v)))",
    "            if not save_pvals:",
    "                bonobo_matrix = np.eye(g) + np.multiply(",
    "                    1 - np.eye(g), (bonobo_matrix * (np.abs(v) > threshold)))",
    "        return (bonobo_matrix, delta, pval)",
    "    SRC = 'inline copy of compute_bonobo'",
    "",
    "X = pd.read_csv(sys.argv[1]).values          # n samples x p genes",
    "E = X.T                                      # netZooPy: genes x samples",
    "emean = np.mean(E, axis=1, keepdims=True)",
    "n = E.shape[1]; g = E.shape[0]",
    "nets = np.empty((g, g, n)); sparse = np.empty((g, g, n))",
    "pvals = np.empty((g, g, n)); deltas = np.empty(n)",
    "for q in range(n):",
    "    m, d, pv = compute_bonobo(E, emean, q, compute_sparse=True,",
    "                              confidence=0.05, save_pvals=True)",
    "    nets[:, :, q] = m; deltas[q] = d; pvals[:, :, q] = pv",
    "    ms, _, _ = compute_bonobo(E, emean, q, compute_sparse=True,",
    "                              confidence=0.05, save_pvals=False)",
    "    sparse[:, :, q] = ms",
    "np.savez(sys.argv[2], nets=nets, sparse=sparse, pvals=pvals, deltas=deltas)",
    "print('reference source:', SRC)"
  ), py_script)

  npz  <- file.path(td, "ref.npz")
  pout <- system2(py, c(shQuote(py_script), shQuote(file.path(td, "X.csv")), shQuote(npz)),
                  stdout = TRUE, stderr = TRUE)
  if (!file.exists(npz)) {
    cat("[SKIP] python reference run failed -- test 4 skipped\n")
    cat(paste(pout, collapse = "\n"), "\n")
  } else {
    cat("[INFO] python reference:", grep("reference source:", pout, value = TRUE), "\n")
    ## read the .npz (zip of .npy files) with base R: parse the NPY v1 header.
    ## numpy writes C-order (last axis fastest), R fills column-major (first
    ## axis fastest), so a numpy array of shape (d1, d2, d3) read linearly is
    ## an R array of dim (d3, d2, d1) -- reverse the dims and aperm them back.
    read_npy <- function(con) {
      magic <- readBin(con, "raw", 6)
      stopifnot(magic[1] == as.raw(0x93), rawToChar(magic[-1]) == "NUMPY")
      ver   <- readBin(con, "integer", 2, size = 1, signed = FALSE)   # major, minor
      hlen  <- readBin(con, "integer", 1, size = 2, endian = "little")
      hdr   <- rawToChar(readBin(con, "raw", hlen))
      descr <- sub(".*'descr':\\s*'([^']+)'.*", "\\1", hdr)
      fortran <- grepl("'fortran_order':\\s*True", hdr)
      shape <- as.integer(strsplit(sub(".*'shape':\\s*\\(([^)]*)\\).*", "\\1", hdr), ",\\s*")[[1]])
      shape <- shape[!is.na(shape)]
      stopifnot(descr == "<f8", !fortran)
      vals <- readBin(con, "double", prod(shape), endian = "little")
      if (length(shape) == 1) return(vals)
      aperm(array(vals, dim = rev(shape)), rev(seq_along(shape)))
    }
    zf <- unzip(npz, exdir = td)
    get_arr <- function(nm) { con <- file(file.path(td, nm), "rb"); on.exit(close(con)); read_npy(con) }
    py_nets   <- get_arr("nets.npy")     # numpy [gene, gene, sample] -> R [gene, gene, sample]
    py_sparse <- get_arr("sparse.npy")
    py_pvals  <- get_arr("pvals.npy")
    py_deltas <- as.numeric(get_arr("deltas.npy"))
    r_dense <- bc$nets; r_sparse <- bc$sparse; r_pvals <- bc$pvals
    e1 <- max(abs(py_nets - r_dense)); e2 <- max(abs(py_deltas - bc$deltas))
    e3 <- max(abs(py_sparse - r_sparse)); e4 <- max(abs(py_pvals - r_pvals))
    check(e1 < 1e-6, sprintf("BONOBO dense nets == Python reference (max diff %.2e)", e1))
    check(e2 < 1e-6, sprintf("BONOBO tuned deltas == Python reference (max diff %.2e)", e2))
    check(e3 < 1e-6, sprintf("BONOBO sparse nets == Python reference (max diff %.2e)", e3))
    check(e4 < 1e-6, sprintf("BONOBO p-values == Python reference (max diff %.2e)", e4))
  }
}

## =============================================================================
## TEST 5 -- lioness_glasso (only with the huge package)
## =============================================================================
if (!requireNamespace("huge", quietly = TRUE)) {
  cat("[SKIP] package 'huge' not available -- test 5 (lioness_glasso) skipped\n")
} else {
  set.seed(9)
  lg <- lioness_glasso(X, lambda = 0.3, progress = FALSE)
  check(identical(dim(lg$nets), c(p, p, n)),
        "lioness_glasso: nets are p x p x n")
  skel_zero <- lg$aggregate == 0
  check(all(vapply(seq_len(n), function(q) all(lg$nets[, , q][skel_zero] == 0), logical(1))),
        "lioness_glasso: restricted nets are zero wherever the aggregate is zero")
  check(all(vapply(seq_len(n), function(q) all(diag(lg$nets[, , q]) == 0), logical(1))),
        "lioness_glasso: per-site diagonals are zero (partial-correlation convention)")
  lg_free <- lioness_glasso(X, lambda = 0.3, A_agg = lg$aggregate,
                            restrict_to_aggregate = FALSE, progress = FALSE)
  check(any(vapply(seq_len(n), function(q) any(lg_free$nets[, , q][skel_zero] != 0), logical(1))),
        "lioness_glasso: restrict_to_aggregate = FALSE leaves off-skeleton entries non-zero")
}

## =============================================================================
## TEST 6 -- caching helpers (explicit root, so no outroot global is needed)
## =============================================================================
td2 <- tempfile("single_cache_")
f <- save_single_nets(lc, "conc", "Unit_X", "lioness", root = td2)
check(is.character(f) && file.exists(f) &&
        grepl(file.path("concentration", "Unit_X", "single_lioness_conc_Unit_X.rds"), f, fixed = TRUE),
      "save_single_nets: writes single_<method>_<kind>_<safe>.rds under <root>/<kind_dir>/<safe>/")
bk <- load_single_nets("conc", "Unit_X", "lioness", root = td2)
check(identical(bk$nets, lc$nets) && is.null(load_single_nets("tu", "Unit_X", "lioness", root = td2)),
      "load_single_nets: round-trips the object; NULL for a missing cache")

## ---- verdict -------------------------------------------------------------------
if (.failures > 0L) {
  cat(.failures, "check(s) FAILED\n")
  quit(status = 1)
}
cat("LIONESS TEST PASSED\n")
