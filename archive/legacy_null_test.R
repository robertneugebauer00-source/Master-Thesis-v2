## ARCHIVED 14.09.2026 -- the legacy Erdos-Renyi-only null.
## Never called: null_test_multi() in scripts/null_models_vRN.R superseded it on
## 27.08.2026. Kept so pre-27.08.2026 numbers stay reproducible. Do NOT source
## this in new work.

## [vRN 27.08.2026] SUPERSEDED for Table 3 by null_test_multi() in
## scripts/null_models_vRN.R. Kept unchanged so every number reported before
## 27.08.2026 stays reproducible. Do NOT call it for new results: it uses ER
## only, which does not preserve the degree sequence and therefore tests a
## different H0 than modularity() assumes internally.
## [vPI2] WEIGHTED Erdos-Renyi null: keep n and density, but give the random graph the
## SAME edge-weight distribution as the observed network, so observed vs null modularity
## are computed on the same (weighted) scale. Average path length is left unweighted on
## both sides (weight would need a distance interpretation).
## @CAVEAT  LEGACY ER-only null, no longer called. Kept for reproducibility only ----------------------
null_test <- function(A, membership, N = 100, seed = 42) {
  ## [vRN-fix 12.08.2026] keep the internal seed for reproducibility, but do NOT leak it
  ## into the global RNG. Without this, any sample() called after a fit restarts from the
  ## same state, so loops that resample sites between fits silently reuse one draw.
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) stats::runif(1)
  .old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit(assign(".Random.seed", .old_seed, envir = .GlobalEnv), add = TRUE)
  g  <- graph_from_adjacency_matrix(abs(A), mode = "undirected", weighted = TRUE, diag = FALSE)
  n  <- vcount(g); p <- edge_density(g)
  w_obs <- E(g)$weight

  blank <- data.frame(metric = c("Modularity","Avg path length"),
                      observed = NA, null_mean = NA, null_sd = NA, z = NA, p_emp = NA)
  if (ecount(g) < 3 || n < 4) return(blank)

  mi      <- as.integer(factor(membership[V(g)$name]))
  obs_mod <- modularity(g, mi, weights = w_obs)             ## [vPI2] weighted
  lcc     <- induced_subgraph(g, which(components(g)$membership == which.max(components(g)$csize)))
  obs_apl <- if (vcount(lcc) > 1) mean_distance(lcc, directed = FALSE, weights = NA) else NA_real_

  set.seed(seed)
  rand_mod <- rand_apl <- numeric(N)
  for (i in seq_len(N)) {
    r <- sample_gnp(n, p)
    if (ecount(r) > 0)
      E(r)$weight <- sample(w_obs, ecount(r), replace = TRUE)  ## [vPI2] resampled observed weights
    rand_mod[i] <- if (ecount(r) > 0)
      modularity(r, membership(cluster_fast_greedy(r)), weights = E(r)$weight) else NA_real_
    cc <- components(r)
    rl <- induced_subgraph(r, which(cc$membership == which.max(cc$csize)))
    rand_apl[i] <- if (vcount(rl) > 1) mean_distance(rl, directed = FALSE, weights = NA) else NA_real_
  }

  z     <- function(obs, x) (obs - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  p_emp <- function(obs, x) {
    x <- x[!is.na(x)]
    (1 + sum(abs(x - mean(x)) >= abs(obs - mean(x)))) / (length(x) + 1)
  }

  data.frame(metric    = c("Modularity", "Avg path length"),
             observed  = c(obs_mod, obs_apl),
             null_mean = c(mean(rand_mod, na.rm = TRUE), mean(rand_apl, na.rm = TRUE)),
             null_sd   = c(sd(rand_mod, na.rm = TRUE),   sd(rand_apl, na.rm = TRUE)),
             z         = c(z(obs_mod, rand_mod), z(obs_apl, rand_apl)),
             p_emp     = c(p_emp(obs_mod, rand_mod), p_emp(obs_apl, rand_apl)))
}
