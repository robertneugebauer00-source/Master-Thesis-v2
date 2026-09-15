## =============================================================================
## null_models_vRN.R  --  modularity null models for the river chemical networks
## Robert Neugebauer, master thesis. Created 19.08.2026.
##
## Replaces the single Erdos-Renyi null in `null_test()` with an explicit,
## switchable hierarchy of null hypotheses. Source this AFTER the workflow's
## function chunks (it needs nothing from them except igraph + huge).
##
##   source("scripts/null_models_vRN.R")
##
## The four nulls, in increasing order of how much of the observed structure
## they hold fixed:
##
##   "er"          Erdos-Renyi G(n, p): fixes only n and global density.
##                 == the null currently used in Clean_Workflow_vRN.Rmd.
##                 H0: "the graph is a random graph of the same size/density."
##                 WARNING: does not preserve the degree sequence. On graphs
##                 with many isolates and heavy degree heterogeneity it is
##                 biased AGAINST finding modularity (see Rhine).
##
##   "config"      Configuration model via degree-preserving rewiring.
##                 H0: "modularity is no more than the degree sequence implies."
##                 This is the null that Newman's Q is *defined* against, so it
##                 is the internally consistent choice. Recommended default for
##                 a graph-level test.
##
##   "weightperm"  Topology fixed, edge weights reshuffled.
##                 H0: "the partial-correlation magnitudes carry no module
##                 information beyond the binary topology."
##
##   ... plus a whole FAMILY of generative null topologies, via gen_null() /
##   null_family_test(): "gnp", "gnm", "config", "chunglu", "ba", "ws10",
##   "ws01", "ring". See the block near the bottom of this file. Short version:
##   degree-aware nulls (config, chunglu) are the defensible ones; lattice-based
##   nulls (ring, ws01, ws10) are built with Q ~ 0.55-0.79 and will call any
##   empirical network "not modular", so do not use them as the headline test.
##
##   "pipeline"    Data-level null: permute the concentration matrix and re-run
##                 the ENTIRE estimator (glasso + StARS + clustering). This is
##                 the only null that accounts for structure manufactured by
##                 the estimator itself. Two flavours, see pipeline_null():
##                   type = "indep" : every compound independent
##                   type = "grad"  : rank-k latent gradient KEPT, residual
##                                    co-structure destroyed  <- the informative one
##
## [repo 15.09.2026] the generative families (gen_null / null_family_test /
## plot_null_grid) and the pipeline null are now WIRED IN: they run in
## analysis/S4_null_hierarchy.Rmd. Only null_test_multi() is used by the
## mainline; everything else is S4 material.
##
## IMPORTANT: Q falls steeply with density, so a pipeline null that returns
## far fewer edges than the observed network cannot be compared on raw Q.
## Use density_matched_Q() / the `match_density` argument.
## =============================================================================

suppressPackageStartupMessages({
  library(igraph)
  library(huge)
})

## ---------------------------------------------------------------------------
## helpers
## ---------------------------------------------------------------------------

.g_from_A <- function(A) {
  graph_from_adjacency_matrix(abs(A), mode = "undirected", weighted = TRUE, diag = FALSE)
}

## same partition the workflow uses; defined here only if the Rmd's version is
## not already in scope, so this file can be sourced standalone by the backfill
## script. Weights passed explicitly -- never rely on the igraph default.
if (!exists("greedy_modules")) {
  greedy_modules <- function(A) {
    g <- .g_from_A(A)
    setNames(as.integer(membership(cluster_fast_greedy(g, weights = E(g)$weight))), V(g)$name)
  }
}

## weighted modularity of the fast-greedy partition (explicit weights: igraph's
## default changed across versions, never rely on weights = NULL)
Q_weighted <- function(g, memb = NULL) {
  if (ecount(g) < 3) return(NA_real_)
  w  <- E(g)$weight
  mi <- if (is.null(memb)) as.integer(membership(cluster_fast_greedy(g, weights = w)))
        else               as.integer(factor(memb[V(g)$name]))
  modularity(g, mi, weights = w)
}

Q_unweighted <- function(g, memb = NULL) {
  if (ecount(g) < 3) return(NA_real_)
  mi <- if (is.null(memb)) as.integer(membership(cluster_fast_greedy(g, weights = NA)))
        else               as.integer(factor(memb[V(g)$name]))
  modularity(g, mi, weights = NA)
}

## keep the m strongest edges, then re-cluster. Use to compare networks of
## different density on equal terms.
thin_graph <- function(g, m_target) {
  if (ecount(g) <= m_target) return(g)
  keep <- order(abs(E(g)$weight), decreasing = TRUE)[seq_len(m_target)]
  subgraph.edges(g, keep, delete.vertices = FALSE)
}

density_matched_Q <- function(A, m_target) {
  g <- thin_graph(.g_from_A(A), m_target)
  c(edges = ecount(g), Q = Q_weighted(g))
}

## ---------------------------------------------------------------------------
## graph-level null test
## ---------------------------------------------------------------------------
## A          weighted adjacency (partial correlations), as returned by fit_glasso()
## membership named integer vector, as returned by greedy_modules()
## nullModel  "config" (default) | "er" | "weightperm"
## weighted   TRUE  -> weighted Q on both sides (matches the vRN convention)
##            FALSE -> unweighted Q on both sides
##
## Returns one row per metric, with the null model recorded in the table so the
## choice is visible in every output file.
## nullModel: "config" | "chunglu" | "er" (= gnp) | "gnm" | "ba" | "kreg" |
##            "grg" | "ws10" | "ws01" | "ring" | "weightperm"
## weightMode: "resample" -> null edges get weights resampled from the observed
##               weight vector. Use this whenever you compare NULL FAMILIES to
##               each other, so that only the topology model differs.
##             "carry"    -> weights travel with the edges (config/weightperm
##               only). Preserves the weight-degree relationship; gives a larger
##               z. Both are defensible; report which you used.
null_test2 <- function(A, membership, N = 200, seed = 42,
                       nullModel = c("config", "chunglu", "er", "gnm", "ba", "kreg",
                                     "grg", "ws10", "ws01", "ring", "weightperm"),
                       weighted   = TRUE,
                       weightMode = c("resample", "carry"),
                       rewire_factor = 20) {

  nullModel  <- match.arg(nullModel)
  weightMode <- match.arg(weightMode)

  ## do not leak the internal seed into the global RNG  [vRN-fix 12.08.2026]
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) stats::runif(1)
  .old <- get(".Random.seed", envir = .GlobalEnv)
  on.exit(assign(".Random.seed", .old, envir = .GlobalEnv), add = TRUE)

  g <- .g_from_A(A)
  n <- vcount(g); m <- ecount(g); p <- edge_density(g); w <- E(g)$weight

  blank <- data.frame(nullModel = nullModel,
                      metric = c("Modularity", "Avg path length"),
                      observed = NA, null_mean = NA, null_sd = NA, z = NA, p_emp = NA)
  if (m < 3 || n < 4) return(blank)

  Qfun    <- if (weighted) Q_weighted else Q_unweighted
  obs_mod <- Qfun(g, membership)

  lcc     <- induced_subgraph(g, which(components(g)$membership == which.max(components(g)$csize)))
  obs_apl <- if (vcount(lcc) > 1) mean_distance(lcc, directed = FALSE, weights = NA) else NA_real_

  set.seed(seed)
  rand_mod <- rand_apl <- numeric(N)

  for (i in seq_len(N)) {
    r <- if (nullModel == "weightperm") g else gen_null(g, if (nullModel == "er") "gnp" else nullModel)
    ## weights: resample by default so families are comparable; "carry" keeps them
    ## attached to the rewired edges (only meaningful for config / weightperm).
    if (ecount(r) > 0) {
      if (nullModel == "weightperm")      E(r)$weight <- sample(w)
      else if (weightMode == "resample")  E(r)$weight <- sample(w, ecount(r), replace = TRUE)
      else if (is.null(E(r)$weight))      E(r)$weight <- sample(w, ecount(r), replace = TRUE)
    }

    rand_mod[i] <- Qfun(r)                    # partition re-optimised on the null graph
    cc <- components(r)
    rl <- induced_subgraph(r, which(cc$membership == which.max(cc$csize)))
    rand_apl[i] <- if (vcount(rl) > 1) mean_distance(rl, directed = FALSE, weights = NA) else NA_real_
  }

  zf <- function(o, x) (o - mean(x, na.rm = TRUE)) / stats::sd(x, na.rm = TRUE)
  pf <- function(o, x) { x <- x[!is.na(x)]
                         (1 + sum(abs(x - mean(x)) >= abs(o - mean(x)))) / (length(x) + 1) }

  data.frame(nullModel = nullModel,
             metric    = c("Modularity", "Avg path length"),
             observed  = c(obs_mod, obs_apl),
             null_mean = c(mean(rand_mod, na.rm = TRUE), mean(rand_apl, na.rm = TRUE)),
             null_sd   = c(stats::sd(rand_mod, na.rm = TRUE), stats::sd(rand_apl, na.rm = TRUE)),
             z         = c(zf(obs_mod, rand_mod), zf(obs_apl, rand_apl)),
             p_emp     = c(pf(obs_mod, rand_mod), pf(obs_apl, rand_apl)))
}

## [repo 15.09.2026] DEPRECATED. Never called by any script in the repo --
## null_test_multi() (config/chunglu/er, degree-preserving first) is the
## drop-in actually used by run_unit(). Kept, with a warning, only so old
## interactive sessions fail loudly instead of silently running the old
## er/config/weightperm set.
null_test_all <- function(A, membership, N = 200, seed = 42, weighted = TRUE) {
  .Deprecated("null_test_multi",
              msg = paste("null_test_all() is deprecated and was never part of the",
                          "reported pipeline -- use null_test_multi() (config/chunglu/er)."))
  do.call(rbind, lapply(c("er", "config", "weightperm"),
    function(nm) null_test2(A, membership, N = N, seed = seed,
                            nullModel = nm, weighted = weighted)))
}

## ---------------------------------------------------------------------------
## DROP-IN for the workflow's null_test()
## ---------------------------------------------------------------------------
## Returns exactly the columns the vRN pipeline already writes to
## table3_nullmodel_sensitivity.csv -- observed, null_mean, null_sd, z, p_emp --
## plus a leading `nullModel` column, one block of rows per null. The rounding
## step in `run_unit()` (num <- c("observed","null_mean","null_sd","z","p_emp"))
## keeps working untouched, so this is a one-word change at the call site:
##
##   cbind(Model = md, null_test(A, memb, N = n_null))
##   ->
##   cbind(Model = md, null_test_multi(A, memb, N = n_null))
##
## Default set, reordered [vRN 05.09.2026]: the two DEGREE-PRESERVING nulls first,
## because they are the H0 that Newman Q is defined against and therefore the headline
## test; ER LAST and supplementary, kept only for continuity with everything reported
## before 31.08.2026 and as the worked counter-example. Values are unchanged -- this
## reorders the rows of table3_nullmodel_sensitivity.csv so the defensible null is read
## first. `weightMode = "resample"` so the three stay comparable to each other.
null_test_multi <- function(A, membership, N = 100, seed = 42,
                            nulls = c("config", "chunglu", "er"),
                            weighted = TRUE, weightMode = "resample") {
  out <- lapply(nulls, function(nm)
    null_test2(A, membership, N = N, seed = seed, nullModel = nm,
               weighted = weighted, weightMode = weightMode))
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}

## ---------------------------------------------------------------------------
## data-level (pipeline) null
## ---------------------------------------------------------------------------
## Xlog   sites x compounds matrix, exactly what you hand to fit_glasso()
## type   "grad"  keep a rank-`rank` latent gradient, destroy the rest
##        "indep" destroy everything
## fitter function(X, seed, rep.num) -> list(A = ...). Defaults to fit_glasso().
##
## Each replicate re-runs glasso + StARS, so cost = R x one model fit.
## ~10 s per fit for a 150-site basin; budget accordingly.
null_data <- function(Xlog, type = c("grad", "indep"), rank = 1) {
  type <- match.arg(type)
  X <- as.matrix(Xlog); dn <- dimnames(X)
  out <- if (type == "indep") {
    apply(X, 2, sample)
  } else {
    ctr <- colMeans(X)
    Z   <- sweep(X, 2, ctr, "-")
    sv  <- svd(Z, nu = rank, nv = rank)
    Fk  <- sv$u[, seq_len(rank), drop = FALSE] %*%
           diag(sv$d[seq_len(rank)], rank, rank) %*%
           t(sv$v[, seq_len(rank), drop = FALSE])
    sweep(Fk + apply(Z - Fk, 2, sample), 2, ctr, "+")
  }
  dimnames(out) <- dn
  out
}

pipeline_null <- function(Xlog, type = c("grad", "indep"), R = 20, rank = 1,
                          rep.num = 20, seed0 = 9000, fitter = NULL, verbose = TRUE) {
  type <- match.arg(type)
  if (is.null(fitter)) {
    stopifnot(exists("fit_glasso"))
    fitter <- function(X, seed, rep.num) fit_glasso(X, seed = seed, rep.num = rep.num)
  }
  do.call(rbind, lapply(seq_len(R), function(i) {
    set.seed(seed0 + i)
    Xn <- null_data(Xlog, type, rank = rank)
    A  <- fitter(Xn, seed = seed0 + i, rep.num = rep.num)$A
    g  <- .g_from_A(A)
    if (verbose) message(sprintf("[null %s] rep %d/%d  E=%d  Q=%.3f",
                                 type, i, R, ecount(g), Q_weighted(g)))
    data.frame(type = type, rep = i, edges = ecount(g), Q = Q_weighted(g))
  }))
}

## Compare an observed network to a pipeline null, on BOTH edge count and
## density-matched modularity. Reporting only raw Q here is misleading.
pipeline_null_test <- function(A_obs, null_df) {
  g_obs  <- .g_from_A(A_obs)
  m_obs  <- ecount(g_obs)
  m_null <- stats::median(null_df$edges)
  m_min  <- min(m_obs, m_null)

  obs_raw     <- Q_weighted(g_obs)
  obs_matched <- Q_weighted(thin_graph(g_obs, m_min))

  zf <- function(o, x) (o - mean(x)) / stats::sd(x)
  data.frame(
    test      = c("Edge count", "Modularity (raw Q)", paste0("Modularity @ ", m_min, " edges")),
    observed  = c(m_obs, obs_raw, obs_matched),
    null_mean = c(mean(null_df$edges), mean(null_df$Q), mean(null_df$Q)),
    null_sd   = c(stats::sd(null_df$edges), stats::sd(null_df$Q), stats::sd(null_df$Q)),
    z         = c(zf(m_obs, null_df$edges), zf(obs_raw, null_df$Q), zf(obs_matched, null_df$Q)),
    note      = c("", "compare only if edge counts are similar", "use this one")
  )
}

## ---------------------------------------------------------------------------
## generative null topologies -- the "which random graph model?" question
## ---------------------------------------------------------------------------
## All families are matched on n and (as closely as the generator allows) on the
## edge count m, then given weights resampled from the observed weight vector,
## so the only thing that varies between families is the TOPOLOGY MODEL.
##
##   gnp      Erdos-Renyi G(n, p)      fixes n + density            [current]
##   gnm      Erdos-Renyi G(n, m)      fixes n + exact edge count
##   config   configuration model      fixes the exact degree sequence
##   chunglu  Chung-Lu                 fixes the EXPECTED degree sequence
##   ba       Barabasi-Albert          scale-free, preferential attachment
##   ws10     Watts-Strogatz p = 0.10  small-world
##   ws01     Watts-Strogatz p = 0.01  near-lattice
##   ring     ring lattice (circular)  every node joined to its nei neighbours
##
## The last three are lattice-derived and carry very high built-in modularity
## (a ring partitions trivially into contiguous arcs). They answer "is my
## network more modular than a circle?", which is not a question anyone needs
## answered. Keep them in a supplementary table as a demonstration, not as a test.

.match_m <- function(gr, m, tries = 60) {
  gr <- simplify(gr, remove.multiple = TRUE, remove.loops = TRUE)
  if (ecount(gr) > m) gr <- delete_edges(gr, sample(seq_len(ecount(gr)), ecount(gr) - m))
  k <- 0
  while (ecount(gr) < m && k < tries) {
    need <- m - ecount(gr); k <- k + 1
    cand <- matrix(sample(vcount(gr), 2 * ceiling(need * 1.4), replace = TRUE), ncol = 2)
    cand <- cand[cand[, 1] != cand[, 2], , drop = FALSE]
    gr <- simplify(add_edges(gr, as.vector(t(cand))))
    if (ecount(gr) > m) gr <- delete_edges(gr, sample(seq_len(ecount(gr)), ecount(gr) - m))
  }
  gr
}

## Which generators have a RING baked in?
##   ring, ws01, ws10   -> yes. igraph's sample_smallworld() *starts* from a
##                         circular lattice and rewires it; make_lattice(circular
##                         = TRUE) is the ring itself. Both are round by default.
##   grg                -> no ring, but still spatially embedded: a random
##                         geometric graph in the unit square. High clustering
##                         WITHOUT a circle -- the right control for showing that
##                         the problem is built-in local structure, not roundness.
##   gnp, gnm, kreg, ba,
##   config, chunglu    -> no geometry at all.
NULL_FAMILIES <- c("gnp", "gnm", "config", "chunglu", "ba", "kreg",
                   "grg", "ws10", "ws01", "ring")

gen_null <- function(g, family) {
  n <- vcount(g); m <- ecount(g); d <- degree(g)
  nei <- max(1L, round(m / n))
  r <- switch(family,
    gnp     = sample_gnp(n, edge_density(g)),
    gnm     = sample_gnm(n, m),
    config  = rewire(g, keeping_degseq(loops = FALSE, niter = 20 * m)),
    chunglu = sample_fitness(m, fitness.out = pmax(d, 1e-8)),
    ba      = sample_pa(n, power = 1, m = nei, directed = FALSE),
    kreg    = { k <- max(1L, round(2 * m / n))
                if ((n * k) %% 2 == 1) k <- k + 1L      # n*k must be even
                sample_k_regular(n, k) },
    grg     = sample_grg(n, sqrt(2 * m / (n * (n - 1) * pi))),  # radius matched to m
    ws10    = sample_smallworld(1, n, nei, p = 0.10),
    ws01    = sample_smallworld(1, n, nei, p = 0.01),
    ring    = make_lattice(length = n, dim = 1, nei = nei, circular = TRUE),
    stop("unknown null family: ", family))
  if (family != "config") r <- .match_m(r, m)
  r
}

## Draw the observed network next to one realisation of each null family, all
## with the same layout algorithm and styling. This is the figure that makes the
## argument without any statistics: the ring nulls literally draw as circles.
plot_null_grid <- function(A, membership, basin, file,
                           families = c("config","chunglu","gnm","ba","kreg",
                                        "grg","ws10","ws01","ring"),
                           seed = 7, width = 2400, height = 1080, res = 190) {
  lab <- c(gnp = "Erdos-Renyi G(n,p)", gnm = "Erdos-Renyi G(n,m)",
           config = "Configuration model", chunglu = "Chung-Lu",
           ba = "Barabasi-Albert", kreg = "Random k-regular",
           grg = "Random geometric", ws10 = "Watts-Strogatz p=0.10",
           ws01 = "Watts-Strogatz p=0.01", ring = "Ring lattice")
  pal <- c("#00549F","#006165","#A11035","#E69F00","#57A0D3","#8E7CC3","#57AB27",
           "#F6A800","#CC071E","#0098A1","#BDCD00","#612158","#7A6FAC","#A9A9A9")
  gobs <- .g_from_A(A)
  drawg <- function(g, ttl) {
    set.seed(seed)
    cl <- membership(cluster_fast_greedy(g, weights = if (is_weighted(g)) E(g)$weight else NA))
    plot(g, layout = layout_with_fr(g), vertex.size = 3.2, vertex.label = NA,
         vertex.color = pal[(as.integer(cl) - 1) %% length(pal) + 1],
         vertex.frame.color = "grey30", vertex.frame.width = 0.4,
         edge.color = adjustcolor("grey35", 0.45), edge.width = 0.5, margin = c(0,0,0,0))
    title(ttl, cex.main = 0.95, line = 0.2)
  }
  png(file, width = width, height = height, res = res)
  on.exit(dev.off(), add = TRUE)
  nr <- ceiling((length(families) + 1) / 5)
  par(mfrow = c(nr, 5), mar = c(0.4, 0.4, 2.2, 0.4))
  drawg(gobs, sprintf("OBSERVED %s\nQ = %.2f", basin, Q_weighted(gobs, membership)))
  for (f in families) {
    set.seed(seed)
    r <- gen_null(gobs, f)
    if (ecount(r) > 0) E(r)$weight <- sample(E(gobs)$weight, ecount(r), replace = TRUE)
    drawg(r, sprintf("%s\nQ = %.2f", lab[f], Q_weighted(r)))
  }
  invisible(file)
}

## Q of the observed network against every family, plus two descriptors that
## show HOW UNLIKE the observed network each null family is (this is what
## justifies the choice, and it is the table to put in the appendix).
null_family_test <- function(A, membership, N = 100, seed = 42,
                             families = NULL_FAMILIES) {
  g <- .g_from_A(A); w <- E(g)$weight
  obs <- Q_weighted(g, membership)
  set.seed(seed)
  do.call(rbind, lapply(families, function(fam) {
    q <- tr <- iso <- numeric(N)
    for (i in seq_len(N)) {
      r <- gen_null(g, fam)
      if (ecount(r) > 0) E(r)$weight <- sample(w, ecount(r), replace = TRUE)
      q[i]   <- Q_weighted(r)
      tr[i]  <- transitivity(r, type = "global")
      iso[i] <- sum(degree(r) == 0)
    }
    data.frame(family = fam, observed = obs,
               null_mean = mean(q), null_sd = stats::sd(q),
               z = (obs - mean(q)) / stats::sd(q),
               p_emp = (1 + sum(q >= obs)) / (N + 1),
               null_clust = mean(tr, na.rm = TRUE), null_isolates = mean(iso))
  }))
}

## =============================================================================
## Example
## =============================================================================
if (FALSE) {

  pu <- prep_unit("conc", sites_meta$site_id[sites_meta$river_basin == "Rhine"])
  gl <- fit_glasso(pu$Xlog, rep.num = 50)
  mb <- greedy_modules(gl$A)

  ## 1. graph-level: all three nulls side by side
  print(null_test_all(gl$A, mb, N = 200))

  ## 2. every generative null family (~5 s for a 150-node network, N = 100)
  print(null_family_test(gl$A, mb, N = 100))
  plot_null_grid(gl$A, mb, "Rhine", "fig_null_graphs_rhine.png")

  ## 3. data-level: gradient-preserving pipeline null  (~10 s x R)
  pn <- pipeline_null(pu$Xlog, type = "grad", R = 20)
  print(pipeline_null_test(gl$A, pn))
}
