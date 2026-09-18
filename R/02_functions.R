## =============================================================================
## 02_functions.R -- preprocessing, estimators, metrics, modules, risk, plotting
## Split out of Clean_Workflow_vRN-14.08.26.Rmd on 14.09.2026 -- code is verbatim.
## Original chunks: fn-prep | fn-glasso-mb | fn-wgcna | #METRIC TABLES (working) etc. | fn-riskdriver | #PLOTTING
## =============================================================================

## ---- fn-prep  [original L743-815] ----

## [vRN 03.09.2026] ---- fixed (global) node set ------------------------------
## Reference pool = every site in an analysed basin. Every spatial unit (country
## group, whole basin, section, reach) is a subset of it, so a compound that
## clears PREV_THRESH here is a candidate node in all of them.
REF_SITES <- sites_meta$site_id[sites_meta$river_basin %in% big_basins]

global_nodeset <- function(kind, prev = PREV_THRESH) {
  mat   <- if (kind == "conc") conc_mat else tu_mat
  sites <- intersect(REF_SITES, colnames(mat))
  sub   <- mat[, sites, drop = FALSE]
  ## @KEY  the prevalence filter itself ---------------------------------------------------------------
  rownames(sub)[rowSums(sub > 0) >= max(2L, ceiling(prev * ncol(sub)))]
}
KEEP_GLOBAL <- list(conc = global_nodeset("conc"), tu = global_nodeset("tu"))
message(sprintf("[vRN] global node set at %.0f %% of %d reference sites: conc %d/%d, tu %d/%d",
                100 * PREV_THRESH, length(intersect(REF_SITES, colnames(conc_mat))),
                length(KEEP_GLOBAL$conc), nrow(conc_mat),
                length(KEEP_GLOBAL$tu),   nrow(tu_mat)))

prep_unit <- function(kind, site_ids, prev = PREV_THRESH, scope = NODESET_SCOPE) {
  mat   <- if (kind == "conc") conc_mat else tu_mat
  sites <- intersect(site_ids, colnames(mat))
  sub   <- mat[, sites, drop = FALSE]

  ## [vRN 03.09.2026] node selection -------------------------------------------
  ## scope = "global": fixed pool (KEEP_GLOBAL), identical for every unit, so
  ##   Q / edge count / density are compared over the same candidate compounds.
  ##   The only within-unit rule left is a zero-variance guard: a compound with
  ##   fewer than MIN_DETECT_UNIT detections here is a constant column once the
  ##   non-detects are floored, and cor() would return NA for it.
  ## scope = "unit": pre-03.09.2026 behaviour, PREV_THRESH of THIS unit sites.
  n_det <- rowSums(sub > 0)
  if (identical(scope, "global")) {
    in_pool <- rownames(sub) %in% KEEP_GLOBAL[[kind]]
    keep    <- in_pool & n_det >= MIN_DETECT_UNIT
    n_guard <- sum(in_pool) - sum(keep)
  } else {
    keep    <- n_det >= max(2L, ceiling(prev * ncol(sub)))
    n_guard <- NA_integer_
  }
  sub <- sub[keep, , drop = FALSE]
  if (!nrow(sub)) stop("prep_unit: no compounds left after node selection")

  ## [vPI3] non-detect substitution: replace only the true zeros (MEC = 0, i.e. below MDL)
  ## with a compound-specific MDL/2 floor. Genuinely detected values are left exactly as
  ## measured -- unlike a pseudocount added to every cell, this does not inflate real
  ## detections (which matters for TU, since that feeds the risk-driver threshold).
  ## conc: mdl_half_ngL, keyed by compound name (built above from MDL_CECs_summary.csv).
  ## tu:   mdl_half_tu where available (MDL/2 converted via EC10); compounds without a
  ##       matched MDL+EC10 fall back to the vPI2 within-unit 0.5*min(nonzero) floor.
  floor_vec <- if (kind == "conc") {
    mdl_half_ngL[rownames(sub)]
  } else {
    fallback_tu <- 0.5 * min(sub[sub > 0])
    v <- mdl_half_tu[rownames(sub)]
    v[is.na(v)] <- fallback_tu
    v
  }
  is_zero     <- sub == 0
  floor_mat   <- matrix(floor_vec, nrow = nrow(sub), ncol = ncol(sub))   # recycled per compound (row)
  sub_floored <- sub
  sub_floored[is_zero] <- floor_mat[is_zero]

  ## @KEY  the matrix every estimator sees: log10, sites x compounds ----------------------------------
  Xlog <- log10(t(sub_floored))   # rows = sites, cols = compounds

  list(cs = sub, Xlog = Xlog, n_sites = ncol(sub), n_comp = nrow(sub),
       pseudo = floor_vec, scope = scope, n_dropped_guard = n_guard)
}


## ---- fn-glasso-mb  [original L929-978] ----
fit_glasso <- function(Xlog, seed = 42, rep.num = STARS_REPNUM, nlambda = 30) {
  ## [vRN-fix 12.08.2026] keep the internal seed for reproducibility, but do NOT leak it
  ## into the global RNG. Without this, any sample() called after a fit restarts from the
  ## same state, so loops that resample sites between fits silently reuse one draw.
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) stats::runif(1)
  .old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit(assign(".Random.seed", .old_seed, envir = .GlobalEnv), add = TRUE)
  set.seed(seed)
  ## @KEY  glasso path: a whole family of networks at increasing penalty ------------------------------
  h   <- huge(as.matrix(Xlog), method = "glasso", nlambda = nlambda,
              cov.output = TRUE, verbose = FALSE)
  sel <- huge.select(h, criterion = "stars", stars.thresh = 0.05,
                     rep.num = rep.num, verbose = FALSE)

  Omega <- as.matrix(sel$opt.icov)          # precision matrix at the chosen lambda
  d     <- sqrt(diag(Omega))
  pcor  <- -Omega / outer(d, d)             # partial correlations
  diag(pcor) <- 0

  ## @KEY  THE PRIMARY NETWORK: partial correlations, StARS-selected edges only -----------------------
  A <- pcor * (as.matrix(sel$refit) != 0)   # keep only the edges StARS selected
  dimnames(A) <- list(colnames(Xlog), colnames(Xlog))
  ## [vRN 05.09.2026] stab = fraction of StARS subsamples in which each pair was
  ## selected, at the chosen lambda. Free (huge.select already computes it); use it
  ## to report edge-level stability, not just the graph-level lambda.
  stab <- as.matrix(sel$merge[[sel$opt.index]]); dimnames(stab) <- dimnames(A)
  list(A = (A + t(A)) / 2, lambda = sel$opt.lambda, stab = stab, rep.num = rep.num)
}

fit_mb    <- function(Xlog, seed = 42, rep.num = STARS_REPNUM, nlambda = 30) {
  ## [vRN-fix 12.08.2026] keep the internal seed for reproducibility, but do NOT leak it
  ## into the global RNG. Without this, any sample() called after a fit restarts from the
  ## same state, so loops that resample sites between fits silently reuse one draw.
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) stats::runif(1)
  .old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit(assign(".Random.seed", .old_seed, envir = .GlobalEnv), add = TRUE)
  set.seed(seed)
  h   <- huge(as.matrix(Xlog), method = "mb", nlambda = nlambda, verbose = FALSE)
  sel <- huge.select(h, criterion = "stars", stars.thresh = 0.05,
                     rep.num = rep.num, verbose = FALSE)

  R <- suppressWarnings(cor(as.matrix(Xlog)))   # marginal sign/weight (cross-check only)
  A <- R * (as.matrix(sel$refit) != 0)
  diag(A) <- 0
  dimnames(A) <- list(colnames(Xlog), colnames(Xlog))
  stab <- as.matrix(sel$merge[[sel$opt.index]]); dimnames(stab) <- dimnames(A)
  list(A = (A + t(A)) / 2, lambda = sel$opt.lambda, stab = stab, rep.num = rep.num)
}

## ---- fn-wgcna  [original L982-1027] ----
## [vPI2] signed-network minimum power scaled to sample size (WGCNA authors' guidance)
min_power_signed <- function(n) if (n < 20) 18 else if (n < 30) 16 else if (n < 40) 14 else 12

fit_wgcna <- function(Xlog, seed = 42, edge_top_frac = 0.05, min_module = 5) {
  X      <- as.matrix(Xlog)
  powers <- 1:20

  ## [vRN 05.09.2026] SOFT THRESHOLD -- PROVENANCE (meeting TODO 04.09.2026).
  ## The soft threshold is Zhang & Horvath (2005, Stat Appl Genet Mol Biol 4:17), written
  ## for MICROARRAY GENE EXPRESSION: raise the correlation to a power beta chosen by the
  ## scale-free-topology criterion, so weak correlations are attenuated rather than cut.
  ## It is NOT a prevalence filter and it is NOT related to the glasso sparsification --
  ## it lives entirely inside this WGCNA branch and touches neither fit_glasso() nor
  ## fit_mb(). The meeting slide listed the two together; they are separate objects
  ## (edge weighting here, node selection there) and should be separated in the methods.
  sft     <- pickSoftThreshold(X, powerVector = powers, networkType = "signed", verbose = 0)
  fit_idx <- -sign(sft$fitIndices$slope) * sft$fitIndices$SFT.R.sq

  above <- which(fit_idx >= 0.80)
  beta  <- if (length(above)) powers[min(above)] else sft$powerEstimate
  if (is.na(beta)) beta <- powers[which.max(fit_idx)]
  beta  <- max(min_power_signed(nrow(X)), min(20, beta))   ## [vPI2] floor scales with n (was max(6, .))

  adj  <- adjacency(X, power = beta, type = "signed")
  diss <- 1 - TOMsimilarity(adj, TOMType = "signed", verbose = 0)
  tree <- hclust(as.dist(diss), method = "average")

  labels  <- cutreeDynamic(tree, distM = diss, deepSplit = 2,
                           pamRespectsDendro = FALSE, minClusterSize = min_module)
  modules <- mergeCloseModules(X, labels2colors(labels), cutHeight = 0.25, verbose = 0)$colors
  names(modules) <- colnames(X)

  # top few % of edges kept for a legible graph (NB: WGCNA density is therefore
  # fixed by construction and is NOT compared to glasso/MB density)
  cut <- quantile(adj[upper.tri(adj)], 1 - edge_top_frac)
  A   <- adj * (adj >= cut)
  diag(A) <- 0
  dimnames(A) <- list(colnames(X), colnames(X))

  list(A = (A + t(A)) / 2, modules = modules, beta = beta,
       sft = sft$fitIndices, tree = tree, merged_colors = modules,
       fit_idx = fit_idx, powers = powers,
       adj = adj, dissTOM = diss)
}

## ---- #METRIC TABLES (working) etc.  [original L1031-1182] ----

metrics_from_adj <- function(A, membership, pool_n = NA_integer_) {
  g <- graph_from_adjacency_matrix(abs(A), mode = "undirected", weighted = TRUE, diag = FALSE)

  if (ecount(g) == 0)
    return(list(g = g, nodes = vcount(g), edges = 0L, density = 0, density_conn = NA_real_,
                mean_clust = NA_real_, modularity = NA_real_, nodes_pool = pool_n,
                avg_path = NA_real_, isolates = vcount(g), lcc_share = NA_real_,
                membership = membership))

  comps <- components(g)
  lcc   <- induced_subgraph(g, which(comps$membership == which.max(comps$csize)))
  mi    <- as.integer(factor(membership[V(g)$name]))

  list(g          = g,
       nodes      = vcount(g),
       edges      = ecount(g),
       density    = edge_density(g),
       density_conn = { keep <- which(degree(g) > 0)                ## [vRN 31.08.2026]
                        if (length(keep) > 1) edge_density(induced_subgraph(g, keep)) else NA_real_ },
       mean_clust = transitivity(g, type = "global"),
       ## @KEY  weighted modularity Q -- the statistic the thesis rests on ----------------------------
       modularity = tryCatch(modularity(g, mi, weights = E(g)$weight),
                             error = function(e) NA_real_),   ## [vPI2] explicit weighted Q
       avg_path   = if (vcount(lcc) > 1) mean_distance(lcc, directed = FALSE, weights = NA) else NA_real_,
       isolates   = sum(degree(g) == 0),                        ## [vRN 27.08.2026]
       lcc_share  = vcount(lcc) / vcount(g),                    ## [vRN 27.08.2026]
       ## [vRN 03.09.2026] `density` is computed on the nodes this unit could actually
       ## fit -- a pool compound needs >= MIN_DETECT_UNIT detections here. That count
       ## still tracks unit size (Spearman cor(n_sites, n_nodes) = 0.92 over the 21
       ## units), so DENSITY IS NOT COMPARABLE ACROSS UNITS and should not be read as
       ## if it were. `nodes_pool` records the pool size so the shortfall is visible
       ## (e.g. Basin_Basque fits 129 of 182). Modularity and the null z ARE comparable:
       ## refitting and then adding the missing compounds back as isolated vertices
       ## leaves both bit-identical, because degree-0 vertices drop out of both terms of
       ## Q and are preserved by keeping_degseq() -- outputs_vRN/isolate_padding_test_20260903.csv.
       nodes_pool = pool_n,
       membership = membership)
}

## [vRN 27.08.2026] weights passed explicitly. igraph's `weights = NULL` default has
## changed meaning across versions -- on 2.2.3 modularity() silently returns the
## UNWEIGHTED Q. Verified: identical partition to the old call on all 17 cached
## units (ARI = 1.000, Q identical to 4 dp), so no reported number changes.
## abs(A) discards the sign of the partial correlations. Checked 27.08.2026:
## 55 of 9452 edges are negative (0.6 %); 9 of 17 units have none. Signed
## modularity is therefore unnecessary -- but say so in the methods.
## [vRN 05.09.2026] WHAT A MODULE MEANS HERE (meeting TODO 04.09.2026: why modularity,
## for which data and hypotheses, and what does a community mean for chemical data?).
##   METHOD ORIGIN. Newman (2006, PNAS 103:8577) frames modularity maximisation as a
##     DATA-ANALYSIS technique, not a hypothesis test: Q is a quality function and the
##     partition is found by optimising it. The paper proposes no significance test, so
##     testing Q against an external ensemble is an addition made downstream, not part
##     of the method (see the Table 3 note).
##   ECOLOGICAL READING. In the microbial literature a module is read as EITHER (i) a
##     group of taxa interacting more with each other than with the rest, OR (ii) a group
##     sharing a niche distinct from other groups (Hernandez et al. 2021, ISME J 15:1722;
##     Hou et al. 2024; "clusters have been interpreted as niches", Faust & Raes 2012,
##     Nat Rev Microbiol 10:538). The hypothesis attached to it is a STABILITY one:
##     modular organisation insulates groups from disturbance, so modularity is predicted
##     to FALL along a stress gradient -- the SGH prediction this thesis tests.
##   CHEMICAL TRANSLATION. Compounds do not interact and have no niche, so neither (i)
##     nor (ii) transfers. An edge here is co-variation of concentration across sites
##     after conditioning, so a module is a set of compounds whose concentrations move
##     together across the river network. Candidate generators: shared EMISSION SOURCE
##     (discharge type / use group), shared TRANSPORT AND DEGRADATION behaviour
##     (parent/metabolite chains, persistence), shared APPLICATION PATTERN (crop calendar,
##     catchment land use).
##   WHICH ONE IS IT -- ANSWERED EMPIRICALLY. A30 (scripts/moa_test_vRN_20260831.R) fits
##     an edge-level nested deviance model on ~18 500 compound pairs per basin: SOURCE is
##     overwhelming (deviance 12.1 to 269.2 on 1 df, p 5.0e-04 to 1.7e-60) and MODE OF
##     ACTION adds nothing on top of it (smallest p = 0.033 over nine tests, with a
##     NEGATIVE coefficient). A20/A29 add the parent/metabolite edge enrichment.
##     So modules in this dataset are EMISSION-SOURCE compartments, not mechanistic or
##     ecological ones. That is a positive finding, not a caveat, and it is the sentence
##     the thesis should carry when it interprets Table 2.
## @KEY  the partition: greedy modularity maximisation on abs(A) --------------------------------------
greedy_modules <- function(A) {
  g <- graph_from_adjacency_matrix(abs(A), mode = "undirected", weighted = TRUE, diag = FALSE)
  setNames(as.integer(membership(cluster_fast_greedy(g, weights = E(g)$weight))), V(g)$name)
}

align_modules <- function(all_names, memb) {
  out <- setNames(rep(NA_integer_, length(all_names)), all_names)
  mi  <- as.integer(factor(memb))
  out[names(memb)] <- mi
  gaps <- which(is.na(out))
  if (length(gaps)) out[gaps] <- max(mi, 0) + seq_along(gaps)
  out
}

ari <- function(a, b)
  igraph::compare(as.integer(factor(a)), as.integer(factor(b)), method = "adjusted.rand")


## ---- fn-riskdriver  [original L1195-1268] ----
cent_from_A <- function(A) {
  g <- graph_from_adjacency_matrix(abs(A), mode = "undirected", weighted = TRUE, diag = FALSE)
  if (ecount(g) == 0) return(setNames(rep(0, ncol(A)), colnames(A)))
  ec <- eigen_centrality(g, weights = E(g)$weight)$vector
  setNames(ec, V(g)$name)
}

## @KEY  centrality and RQ kept as two INDEPENDENT node attributes ------------------------------------
risk_driver_table <- function(site_ids, network_A) {
  sites <- intersect(site_ids, colnames(tu_mat))
  sub   <- tu_mat[, sites, drop = FALSE]
  maxRQ <- apply(sub, 1, max, na.rm = TRUE)
  Nsite <- rowSums(sub >= RQ_THRESH, na.rm = TRUE)
  driver <- maxRQ >= RQ_THRESH

  ## [vRN] share of the unit's TOTAL toxic pressure (sum-TU) carried by each compound,
  ## i.e. "% of the calculated risk" in Pedro's sense (complements the per-site RQ rule).
  sumTU_comp <- rowSums(sub, na.rm = TRUE)
  totTU      <- sum(sumTU_comp, na.rm = TRUE)
  pct_sumTU  <- if (totTU > 0) 100 * sumTU_comp / totTU
                else setNames(rep(0, length(sumTU_comp)), names(sumTU_comp))

  netc <- colnames(network_A)
  cent <- cent_from_A(network_A)
  hub_thr <- if (length(cent) && any(cent > 0)) quantile(cent, 0.90, na.rm = TRUE) else Inf

  # union of network compounds and any risk driver (so filtered-out drivers stay visible)
  keep <- union(netc, names(driver)[driver])
  data.frame(
    Compound      = keep,
    Use_group     = ug_lookup[keep],
    in_network    = keep %in% netc,
    centrality    = round(ifelse(keep %in% names(cent), cent[keep], NA_real_), 4),
    structural_hub = ifelse(keep %in% names(cent), cent[keep] >= hub_thr, FALSE),
    maxRQ         = round(ifelse(keep %in% names(maxRQ), maxRQ[keep], NA_real_), 4),
    pct_sumTU     = round(ifelse(keep %in% names(pct_sumTU), pct_sumTU[keep], NA_real_), 2),  ## [vRN]
    N_driver_sites = ifelse(keep %in% names(Nsite), Nsite[keep], 0L),
    absolute_risk_driver = ifelse(keep %in% names(driver), driver[keep], FALSE),
    row.names = NULL
  )[order(-(keep %in% names(driver) & driver[match(keep, names(driver))])), ]
}

plot_hub_vs_risk <- function(rd, file, unit_label) {
  png(file, width = 1500, height = 1150, res = 150); par(mar = c(4.6, 4.8, 4, 1))
  x <- rd$centrality; x[is.na(x)] <- 0
  y <- pmax(rd$maxRQ, 1e-4); y[is.na(y)] <- 1e-4
  ug <- rd$Use_group; ug[is.na(ug) | !(ug %in% names(ug_palette))] <- "other"
  hub_thr <- suppressWarnings(quantile(rd$centrality[rd$in_network], 0.90, na.rm = TRUE))
  ## [vRN 14.08.2026] drivers outlined in RISK_RING so fig4 speaks the same visual language
  ## as the network panels (violet = absolute risk driver, everywhere).
  plot(x, y, log = "y",
       pch = ifelse(rd$absolute_risk_driver, 24, 21),
       bg  = ug_palette[ug],
       col = ifelse(rd$absolute_risk_driver, RISK_RING, "grey55"),
       lwd = ifelse(rd$absolute_risk_driver, 2.2, 0.8),
       cex = ifelse(rd$absolute_risk_driver, 1.7, 1.0),
       xlab = "Structural centrality (eigenvector, glasso network)",
       ylab = "max RQ across unit sites (TU)",
       main = sprintf("Sentinel vs. absolute risk driver — %s", unit_label))
  abline(h = RQ_THRESH, col = "red", lty = 2)
  if (is.finite(hub_thr)) abline(v = hub_thr, col = "grey50", lty = 3)
  text(par("usr")[1], RQ_THRESH, "RQ = 0.02", col = "red", pos = 4, cex = 0.75, offset = 0.2)
  drv <- which(rd$absolute_risk_driver)
  if (length(drv)) text(x[drv], y[drv], labels = substr(rd$Compound[drv], 1, 20),
                        pos = 3, cex = 0.6, col = "grey15", xpd = NA)
  legend("bottomright", bty = "n", cex = 0.8,
         legend = c("absolute risk driver (violet-outlined triangle)", "other compound (circle)",
                    "risk drivers left of the dotted line are peripheral / off-network"),
         pch = c(24, 21, NA), pt.bg = c("grey70", "grey70", NA),
         col = c(RISK_RING, "grey55", NA), pt.lwd = c(2.2, 0.8, NA))
  dev.off()
}

## ---- #METRIC ROW BUILDER (canonical, single source of truth) ---------------
## [repo 15.09.2026] `row()` used to be defined twice -- in run_unit()
## (03_mainline.Rmd) and in refresh_unit_figures() (99_utils.R) -- and the two
## copies had ALREADY drifted (the 99_utils copy lacked `Nodes (pool)`).
## metric_row() is now the single definition; both call sites use it.
metric_row <- function(m) c(Nodes = m$nodes, Edges = m$edges, Density = round(m$density, 4),
                            ## `Density` counts the degree-0 compounds in n, so it is
                            ## diluted by them -- and it is exactly the p the ER null is
                            ## matched on. `Density (excl. isolates)` is the density of
                            ## the part of the network that actually has edges.
                            `Density (excl. isolates)` = round(m$density_conn, 4),
                            ## `Nodes` is what this unit could fit; `Nodes (pool)` is the
                            ## global pool it was drawn from. Report both; do NOT compare
                            ## `Density` across units -- the denominators differ.
                            `Nodes (pool)` = m$nodes_pool,
                            `Mean clustering (unweighted)` = round(m$mean_clust, 3),
                            `Modularity Q` = round(m$modularity, 3),
                            ## avg_path is measured on the LARGEST COMPONENT only.
                            `Avg path length` = round(m$avg_path, 3),
                            Isolates = m$isolates,
                            `LCC share` = round(m$lcc_share, 3))

## ---- #PLOTTING  [original L1272-1478] ----

mod_palette <- function(memb) {
  base <- c("#E63946","#457B9D","#2A9D8F","#E9C46A","#9B2226","#F4A261","#6A994E","#7678ED",
            "#C77DFF","#48CAE4","#B5838D","#264653","#BC6C25","#A5A58D","#606c38","#bc4749",
            "#3a86ff","#ff006e","#8338ec","#fb5607","#ffbe0b","#118ab2")
  ids  <- sort(unique(memb))
  cols <- if (length(ids) <= length(base)) base[seq_along(ids)] else grDevices::rainbow(length(ids))
  setNames(cols, ids)
}

node_cols_module <- function(A, memb) {
  v   <- memb[colnames(A)]
  pal <- mod_palette(v)
  setNames(pal[as.character(v)], colnames(A))
}

node_cols_ug <- function(A) {
  ug <- ug_lookup[colnames(A)]
  ug[is.na(ug) | !(ug %in% names(ug_palette))] <- "other"
  setNames(ug_palette[ug], colnames(A))
}

ug_present_in <- function(As) {
  seen <- unique(unlist(lapply(As, function(A) {
    ug <- ug_lookup[colnames(A)]
    ug[is.na(ug) | !(ug %in% names(ug_palette))] <- "other"
    ug
  })))
  intersect(names(ug_palette), seen)
}

## [vRN 14.08.2026] risk-driver marker redesigned. The old version drew a red frame directly
## on the node; since ~half of all drivers are pharmaceuticals (fill #E63946) that was
## red-on-red and effectively invisible. The marker is now a coloured ring around the node,
## drawn BEHIND it and flush against its edge, so it reads on any fill and stays
## distinct from the black hub ring. RISK_RING is deliberately a hue no use group uses.
RISK_RING <- "#7B2CBF"   ## violet — furthest in RGB from every colour in ug_palette

draw_net <- function(A, node_cols, title, seed = 42, risk = NULL) {
  g <- graph_from_adjacency_matrix(abs(A), mode = "undirected", weighted = TRUE, diag = FALSE)
  if (ecount(g) == 0) { plot.new(); title(main = title); text(0.5, 0.5, "no stable edges", col = "grey50"); return(invisible()) }

  eig   <- eigen_centrality(g, weights = E(g)$weight)$vector
  hub   <- eig >= quantile(eig, 0.90, na.rm = TRUE)
  is_rd <- if (!is.null(risk)) V(g)$name %in% risk else rep(FALSE, vcount(g))
  set.seed(seed)
  lay <- layout_with_fr(g)            # computed once and reused by all three layers, so the
  vs  <- 3 + 8 * (eig / max(eig))     # node positions are identical to the previous figures
  ew  <- 0.4 * (E(g)$weight / max(E(g)$weight)) + 0.2

  ## layer 1 (back): edges, title, and the driver ring. The ring is only slightly larger
  ## than the node, so it sits flush against the node edge like the hub frame does --
  ## a detached ring reads well but eats far too much space in the dense core.
  ## pmax(.., vs + 2) keeps the band visible on the smallest (least central) nodes.
  plot(g, layout = lay, vertex.size = ifelse(is_rd, pmax(vs * 1.40, vs + 2), 0),
       vertex.color = RISK_RING, vertex.frame.color = NA, vertex.label = NA,
       edge.color = "grey80", edge.width = ew, main = title)
  ## layer 2 (front): the nodes themselves; hubs keep their black frame, so a compound that
  ## is both a sentinel and a driver shows both marks (violet band outside, black frame on).
  plot(g, layout = lay, add = TRUE, vertex.color = node_cols[V(g)$name],
       vertex.frame.color = ifelse(hub, "black", NA),
       vertex.frame.width = ifelse(hub, 2.5, 0),
       vertex.size = vs, vertex.label = NA, edge.color = NA)
}

module_legend <- function(A, memb, method, max_mods = 6) {
  nodes <- colnames(A)
  m     <- memb[nodes]
  ug    <- ug_lookup[nodes]; ug[is.na(ug) | !(ug %in% names(ug_palette))] <- "other"
  pal   <- mod_palette(m)

  sizes <- sort(table(m), decreasing = TRUE)
  ids   <- names(sizes)[seq_len(min(max_mods, length(sizes)))]
  labs  <- sapply(ids, function(id) {
    top <- sort(table(ug[m == id]), decreasing = TRUE)
    sprintf("%d compounds  ·  %s", sizes[[id]],
            paste(names(top)[seq_len(min(2, length(top)))], collapse = ", "))
  })

  plot.new()
  legend("topleft", bty = "n", cex = 0.95, title = sprintf("%s modules", method),
         title.adj = 0, legend = labs, pch = 22, pt.bg = pal[ids],
         pt.cex = 2.1, col = "grey40", y.intersp = 1.15)
  if (length(sizes) > max_mods)
    mtext(sprintf("(+%d smaller modules)", length(sizes) - max_mods),
          side = 1, line = -1, cex = 0.75, col = "grey55")
}

panel_title <- function(metrics, md)
  sprintf("%s\nQ=%.2f | E=%d | dens=%.3f",
          c(Glasso = "Glasso (primary)", MB = "MB (cross-check)", WGCNA = "WGCNA (modules)")[md],
          metrics[[md]]["Modularity Q"], metrics[[md]]["Edges"], metrics[[md]]["Density"])

plot_modules <- function(As, cols, mods, metrics, file, unit_label, kind_label, risk = NULL) {
  png(file, width = 2400, height = 1500, res = 150)
  layout(matrix(c(1,2,3, 4,5,6), nrow = 2, byrow = TRUE), heights = c(5, 2.4))
  par(mar = c(1.5, 1, 4, 1), oma = c(0, 0, 3, 0))
  draw_net(As$Glasso, cols$Glasso, panel_title(metrics, "Glasso"), risk = risk)
  draw_net(As$MB,     cols$MB,     panel_title(metrics, "MB"))
  draw_net(As$WGCNA,  cols$WGCNA,  panel_title(metrics, "WGCNA"))
  par(mar = c(1, 1, 1, 1))
  module_legend(As$Glasso, mods$Glasso, "Glasso")
  module_legend(As$MB,     mods$MB,     "MB")
  module_legend(As$WGCNA,  mods$WGCNA,  "WGCNA")
  mtext(sprintf("%s — %s   (nodes coloured by module; size = eigenvector centrality; black ring = hub; violet ring = absolute risk driver)",
                unit_label, kind_label), outer = TRUE, cex = 1.1, font = 2)
  dev.off()
}

plot_usegroups <- function(As, cols, metrics, file, unit_label, kind_label, ug_present, risk = NULL) {
  png(file, width = 2400, height = 1180, res = 150)
  layout(matrix(c(1,2,3, 4,4,4), nrow = 2, byrow = TRUE), heights = c(5.3, 1.7))
  par(mar = c(1.5, 1, 4, 1), oma = c(0, 0, 3, 0))
  draw_net(As$Glasso, cols$Glasso, panel_title(metrics, "Glasso"), risk = risk)
  draw_net(As$MB,     cols$MB,     panel_title(metrics, "MB"))
  draw_net(As$WGCNA,  cols$WGCNA,  panel_title(metrics, "WGCNA"))
  par(mar = c(0.5, 1, 0.5, 1)); plot.new()
  ## horiz = TRUE silently clips when many use groups are present (pharmaceutical and
  ## pesticide fall off the left edge) -- wrap into two rows instead
  legend("top", ncol = ceiling(length(ug_present) / 2), bty = "n", cex = 1.0,
         title = "Node colour — compound use group",
         legend = ug_present, fill = ug_palette[ug_present], border = NA)
  legend("bottom", horiz = TRUE, bty = "n", cex = 1.05,
         legend = c("Node size = eigenvector centrality", "Hub = top 10% (black ring)",
                    "Violet ring = absolute risk driver", "Edge = association"),
         pch = c(21, 21, 21, NA), pt.bg = c("grey60", "grey60", "white", NA),
         pt.cex = c(1.6, 1.8, 1.8, NA), pt.lwd = c(1, 2.4, 3, NA),
         col = c("grey40", "black", RISK_RING, NA), lty = c(NA, NA, NA, 1))
  mtext(sprintf("%s — %s   (nodes coloured by use group)", unit_label, kind_label),
        outer = TRUE, cex = 1.1, font = 2)
  dev.off()
}

plot_wgcna_diag <- function(wg, file, unit_label, kind_label) {
  comp    <- names(wg$modules)
  ug      <- ug_lookup[comp]; ug[is.na(ug) | !(ug %in% names(ug_palette))] <- "other"
  ug_cols <- setNames(unname(ug_palette[ug]), comp)

  f1 <- tempfile(fileext = ".png"); f2 <- tempfile(fileext = ".png")
  png(f1, width = 1000, height = 900, res = 150); par(mar = c(4.5, 4.5, 4, 1))
  plot(wg$powers, wg$fit_idx, type = "n",
       xlab = "Soft-threshold power (beta)", ylab = "Scale-free topology fit (signed R^2)",
       main = sprintf("Scale-free topology fit\n%s — %s", unit_label, kind_label),
       ylim = c(min(0, min(wg$fit_idx)), 1))
  text(wg$powers, wg$fit_idx, labels = wg$powers,
       col = ifelse(wg$powers == wg$beta, "red", "black"), cex = 0.9)
  abline(h = 0.80, col = "red", lty = 2); abline(v = wg$beta, col = "grey60", lty = 3)
  legend("bottomright", legend = sprintf("chosen beta = %d", wg$beta), bty = "n", text.col = "red")
  dev.off()

  png(f2, width = 1150, height = 900, res = 150)
  plotDendroAndColors(wg$tree, ug_cols[comp], "Use group", dendroLabels = FALSE,
                      hang = 0.03, addGuide = TRUE, guideHang = 0.05,
                      main = sprintf("Compound dendrogram — %s (%s)", unit_label, kind_label))
  dev.off()

  image_write(image_append(c(image_read(f1), image_read(f2))), file)
  file.remove(f1, f2)
}

plot_wgcna_dashboard <- function(wg, file, unit_label, kind_label) {
  adj <- wg$adj; TOM <- 1 - wg$dissTOM; dissTOM <- wg$dissTOM; tree <- wg$tree; modules <- wg$modules
  comp <- colnames(adj)
  dimnames(TOM) <- dimnames(dissTOM) <- list(comp, comp)

  ug <- ug_lookup[comp]; ug[is.na(ug) | !(ug %in% names(ug_palette))] <- "other"
  ucol    <- setNames(unname(ug_palette[ug]), comp)
  present <- intersect(names(ug_palette), unique(ug))

  f_bar <- tempfile(fileext = ".png")
  m <- table(factor(ug, levels = present),
             factor(modules[comp], levels = names(sort(table(modules), decreasing = TRUE))))
  dom  <- apply(m, 2, function(cc) rownames(m)[which.max(cc)])
  labs <- sprintf("%s (%s)", dom, colnames(m))
  png(f_bar, width = 1700, height = 980, res = 150); par(mar = c(12, 4.5, 4, 11))
  barplot(m, col = ug_palette[rownames(m)], border = "grey40", las = 2, names.arg = labs, cex.names = 0.9,
          ylab = "Module size (N compounds)", main = "(e) Module size & use-group composition")
  legend("topright", inset = c(-0.17, 0), xpd = NA, legend = present,
         fill = ug_palette[present], bty = "n", cex = 0.8, title = "Use group"); dev.off()

  f_tom <- tempfile(fileext = ".png"); plotTOM <- dissTOM^7; diag(plotTOM) <- NA
  png(f_tom, width = 1200, height = 1150, res = 150)
  TOMplot(plotTOM, tree, ucol[comp], main = "(c) TOM heatmap (bars = use group)",
          col = colorRampPalette(c("#FFFFCC", "#FD8D3C", "#800026"))(250)); dev.off()

  set.seed(42); nn <- min(15, length(comp) - 1)
  emb <- uwot::umap(TOM, n_neighbors = nn, min_dist = 0.3, verbose = FALSE); rownames(emb) <- comp
  conn <- rowSums(adj); lab_i <- order(conn, decreasing = TRUE)[1:min(12, length(conn))]
  lab  <- ifelse(nchar(comp) > 16, paste0(substr(comp, 1, 15), "…"), comp)
  f_um <- tempfile(fileext = ".png"); png(f_um, width = 1300, height = 1150, res = 150); par(mar = c(4.2, 4.2, 4, 1))
  plot(emb, col = ucol[comp], pch = 19, cex = 1.1, xlab = "UMAP-x", ylab = "UMAP-y",
       main = "(d) Feature UMAP (coloured by use group)")
  text(emb[lab_i, 1], emb[lab_i, 2], labels = lab[lab_i], cex = 0.62, pos = 3, col = "grey15", xpd = NA)
  legend("topright", legend = present, col = ug_palette[present], pch = 19, pt.cex = 1.1,
         cex = 0.7, bty = "n", title = "Use group"); dev.off()

  top <- image_append(c(image_read(f_tom), image_read(f_um)))
  bar <- image_scale(image_read(f_bar), paste0(image_info(top)$width))
  fin <- image_append(c(top, bar), stack = TRUE)
  ttl <- image_annotate(image_blank(image_info(fin)$width, 70, color = "white"),
           sprintf("WGCNA diagnostics — %s, %s", unit_label, kind_label),
           gravity = "center", size = 34, weight = 700, color = "#222222")
  image_write(image_append(c(ttl, fin), stack = TRUE), file)
  file.remove(f_bar, f_tom, f_um)
}


## =============================================================================
## [repo 15.09.2026] ---- fitted-artefact cache + run manifest ------------------
##
## run_unit() (03_mainline.Rmd) now saves the FITTED OBJECTS -- adjacency
## matrices, memberships, StARS lambdas -- as one .rds per unit per data kind.
## Everything downstream (99_utils figure refresh, S2 SBM comparison, S3 basin
## panel) reads these instead of refitting. This is what makes "figure refresh
## without refitting" literally true, and it removes S2's dependency on the
## external vRN_followups.R cache.
## =============================================================================

## "conc" -> "concentration", "tu" -> "toxic_unit" (the two output trees)
kind_dir <- function(kind) if (kind == "conc") "concentration" else "toxic_unit"

fitted_path <- function(kind, safe, root = outroot)
  file.path(root, kind_dir(kind), safe, sprintf("fitted_%s_%s.rds", kind, safe))

## Returns the cached fit list, or NULL if the unit was never run / cache absent.
load_fitted <- function(kind, safe, root = outroot) {
  f <- fitted_path(kind, safe, root)
  if (file.exists(f)) readRDS(f) else NULL
}

## fit: list(A_glasso, A_mb, A_wgcna, mods_glasso, mods_mb, mods_wgcna,
##           lambda_glasso, lambda_mb, stab_glasso, beta, kind, unit, n_sites,
##           n_comp, rep.num)
save_fitted <- function(fit, kind, safe, root = outroot) {
  f <- fitted_path(kind, safe, root)
  dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
  saveRDS(fit, f)
  invisible(f)
}

## [repo 15.09.2026] one manifest per run: every number in the output tree is
## traceable to the exact parameter set, package versions and git commit that
## produced it. JSON when jsonlite is available, plain text otherwise.
write_run_manifest <- function(root = outroot, extra = list()) {
  grab <- function(nm) if (exists(nm, envir = globalenv())) get(nm, envir = globalenv()) else NULL
  si  <- sessionInfo()
  pk  <- c(si$loadedOnly, si$otherPkgs)
  git <- tryCatch(system2("git", c("-C", ".", "rev-parse", "--short", "HEAD"),
                          stdout = TRUE, stderr = FALSE), error = function(e) NULL)
  manifest <- c(
    list(timestamp   = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
         r_version   = R.version.string,
         git_commit  = if (length(git)) git[[1]] else NA_character_,
         parameters  = list(PREV_THRESH = grab("PREV_THRESH"),
                            NODESET_SCOPE = grab("NODESET_SCOPE"),
                            STARS_REPNUM = grab("STARS_REPNUM"),
                            RQ_THRESH = grab("RQ_THRESH"),
                            BASIN_MIN_N = grab("BASIN_MIN_N"),
                            SECTION_MIN_N = grab("SECTION_MIN_N"),
                            MIN_DETECT_UNIT = grab("MIN_DETECT_UNIT"),
                            DEFAULT_PSEUDO_CONC = grab("DEFAULT_PSEUDO_CONC")),
         packages    = vapply(pk, function(p) as.character(p$Version), character(1))),
    extra)
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    jsonlite::write_json(manifest, file.path(root, "run_manifest.json"),
                         auto_unbox = TRUE, pretty = TRUE, null = "null")
  } else {
    saveRDS(manifest, file.path(root, "run_manifest.rds"))
    writeLines(capture.output(str(manifest)), file.path(root, "run_manifest.txt"))
  }
  invisible(manifest)
}


## =============================================================================
## [repo 15.09.2026] ---- module stability under site resampling ----------------
##
## S3 builds co-membership matrices for the rarefied cross-basin comparison;
## this applies the same idea WITHIN one unit: refit the glasso network on R
## site subsamples and report how stable the REPORTED partition is (mean ARI
## against the reference partition, plus the averaged co-membership matrix).
## Every module-level claim can then carry a stability number.
##
## REFERENCE PARTITION: pass `safe` (the unit's safe name) to compare against
## the cached mainline fit (STARS_REPNUM = 100) -- that is the partition the
## thesis reports. Only without a cache is the reference refit here, at
## `rep.num`. The subsample fits all use `rep.num` (fixed across draws, so the
## draws are comparable to each other); read the result as "how much the
## reported partition moves under site resampling", not as an absolute
## agreement between two StARS depths.
## =============================================================================
bootstrap_module_stability <- function(kind, site_ids, R = 25, frac = 0.8,
                                       seed = 42, rep.num = 20, safe = NULL) {
  fc <- if (!is.null(safe)) tryCatch(load_fitted(kind, safe), error = function(e) NULL) else NULL
  if (!is.null(fc)) {
    A_full <- fc$A_glasso; m_full <- fc$mods_glasso
  } else {
    pu_full <- prep_unit(kind, site_ids)
    A_full  <- fit_glasso(pu_full$Xlog, rep.num = rep.num)$A
    m_full  <- greedy_modules(A_full)
  }

  ## [vRN-fix pattern] draw every subsample BEFORE any of the subsample fits,
  ## so the draws are independent of whatever the fitter does to the RNG.
  set.seed(seed)
  draws <- lapply(seq_len(R), function(i)
    sample(site_ids, max(4L, floor(frac * length(site_ids)))))

  comps <- colnames(A_full)
  acc   <- matrix(0, length(comps), length(comps), dimnames = list(comps, comps))
  cnt   <- acc
  aris  <- rep(NA_real_, R)
  for (r in seq_len(R)) {
    pu <- tryCatch(prep_unit(kind, draws[[r]]), error = function(e) NULL)
    if (is.null(pu) || pu$n_comp < 5) next
    A  <- tryCatch(fit_glasso(pu$Xlog, seed = seed + r, rep.num = rep.num)$A,
                   error = function(e) NULL)
    if (is.null(A)) next
    mm <- greedy_modules(A)
    common <- intersect(names(m_full), names(mm))
    if (length(common) >= 5) aris[r] <- ari(m_full[common], mm[common])
    nm <- intersect(names(mm), comps)
    if (length(nm) >= 2) {
      same <- outer(mm[nm], mm[nm], `==`)
      acc[nm, nm] <- acc[nm, nm] + same; cnt[nm, nm] <- cnt[nm, nm] + 1
    }
  }
  list(mean_ari = mean(aris, na.rm = TRUE), sd_ari = stats::sd(aris, na.rm = TRUE),
       ari = aris, comembership = acc / pmax(cnt, 1),
       R_requested = R, R_done = sum(!is.na(aris)), frac = frac, rep.num = rep.num,
       reference = if (!is.null(fc)) "cached mainline fit" else sprintf("refit at rep.num = %d", rep.num))
}


## =============================================================================
## [repo 15.09.2026] ---- signed-edge check -------------------------------------
##
## abs(A) discards the sign of ~0.6 % of edges (checked 27.08.2026: 55 of 9452;
## 9 of 17 units have none). signed_edge_check() turns "we checked" into a
## number: it reports the negative-edge share and the ARI between the partition
## on |A| and the partition on the positive-edges-only network. ARI ~ 1 means
## dropping signs cannot change any module claim.
## =============================================================================
signed_edge_check <- function(A) {
  ut    <- upper.tri(A) & A != 0
  n_neg <- sum(A[ut] < 0)
  n_edg <- sum(ut)
  A_pos <- pmax(A, 0)
  m_abs <- greedy_modules(A)
  m_pos <- greedy_modules(A_pos)
  common <- intersect(names(m_abs), names(m_pos))
  data.frame(edges = n_edg, negative_edges = n_neg,
             pct_negative = round(100 * n_neg / max(1, n_edg), 2),
             ARI_abs_vs_positive = if (length(common) >= 5) round(ari(m_abs[common], m_pos[common]), 3) else NA_real_,
             row.names = NULL)
}

