## =============================================================================
## 10_lenses_vRN.R -- which LENS explains the modules: source vs structure vs MOA
## Robert Neugebauer (vRN). Created 15.09.2026.
##
## WHY THIS FILE EXISTS.
##   metrics_from_adj()'s header records the A30 result: modules here are
##   EMISSION-SOURCE compartments -- source (use group) overwhelming, mode of
##   action adding nothing on top of it. A30 pitted TWO lenses against each
##   other. This file closes the question against a THIRD, structurally
##   independent lens -- ClassyFire chemical family -- and, crucially, scores all
##   available lenses on the SAME networks and the SAME common node set, so
##   "most fitting" is a fair comparison and not an artefact of one lens
##   labelling more compounds than another.
##
##   THREE SCORES PER LENS, from complementary angles:
##     * attribute modularity Q_attr -- weighted modularity of the network when
##       the LENS ITSELF is imposed as the partition (module-level fit), with a
##       node-label PERMUTATION null (Q_z, Q_p) that is honest about the
##       non-independence of compound pairs.
##     * nominal assortativity -- a parameter-free same-label-edge descriptor.
##     * edge-level nested deviance -- the A30 statistic, generalised: ONE joint
##       logistic model of edge presence over all pairs gives each lens a
##       MARGINAL deviance (alone) and a UNIQUE drop-one deviance (what it adds
##       OVER the others). Because the lenses are correlated (same-source
##       compounds often share structure), marginal and unique diverge -- the
##       UNIQUE column is the one that answers "which lens fits best", and the
##       coefficient SIGN flags a lens that is significant but NEGATIVE (as MOA
##       was in A30: same-label pairs LESS connected -- not a module generator).
##
##   Headline = unique deviance; Q_z corroborates at the module level. Coverage
##   is reported per lens because ClassyFire has gaps; the fair ranking uses the
##   node set every lens can label.
##
## DEPENDS ON: 00_setup.R, 01_data.R (ug_lookup, compounds, compound_inchikey),
##   02_functions.R (load_fitted, prep_unit, fit_glasso, greedy_modules). igraph.
## =============================================================================

LENS_MIN_NODES <- 8L       # skip a unit with fewer co-labelled nodes than this
LENS_PERM_B    <- 999L     # permutations for the modularity null
CF_LEVEL       <- "class"  # ClassyFire level for the structure lens: "class" or
                           # "superclass". Switch to "superclass" if class is too
                           # sparse to form >= 2 populated groups (~500 compounds
                           # collapse into a few benzenoid/heterocycle classes).

## ---------------------------------------------------------------------------
## 1. assemble the lenses -- source (ug_lookup), structure (ClassyFire), MOA
## ---------------------------------------------------------------------------
## structure: in-memory cf_lookup (from the S1 {#classyfire} chunk) if present at
## class level, else rebuilt from outroot/classyfire_lookup.csv at CF_LEVEL.
build_cf_lookup <- function(level = CF_LEVEL) {
  if (level == "class" && exists("cf_lookup", inherits = TRUE)) return(get("cf_lookup", inherits = TRUE))
  cache <- file.path(outroot, "classyfire_lookup.csv")
  if (!file.exists(cache)) return(NULL)
  cf <- utils::read.csv(cache, stringsAsFactors = FALSE)
  if (!all(c("InChIKey", level) %in% names(cf))) return(NULL)
  setNames(cf[[level]][match(compound_inchikey, cf$InChIKey)], compounds)
}
## MOA: optional raw_data/compound_moa.csv (compound|inchikey + moa/mechoa/verhaar)
build_moa_lookup <- function() {
  f <- if (exists("in_file")) in_file("compound_moa.csv", "raw_data") else file.path(proj, "raw_data", "compound_moa.csv")
  if (!file.exists(f)) return(NULL)
  m <- utils::read.csv(f, stringsAsFactors = FALSE)
  moa_col <- intersect(c("moa", "MOA", "mode_of_action", "mechoa", "verhaar"), names(m))[1]
  if (is.na(moa_col)) return(NULL)
  ik_col   <- intersect(c("inchikey", "InChIKey"), names(m))[1]
  name_col <- intersect(c("compound", "chemical_name", "name"), names(m))[1]
  out <- setNames(rep(NA_character_, length(compounds)), compounds)
  if (!is.na(ik_col))   out[] <- setNames(m[[moa_col]], toupper(trimws(m[[ik_col]])))[toupper(trimws(compound_inchikey))]
  if (!is.na(name_col)) { miss <- is.na(out)
    out[miss] <- setNames(m[[moa_col]], tolower(trimws(m[[name_col]])))[tolower(trimws(compounds[miss]))] }
  if (all(is.na(out))) NULL else out
}
assemble_lenses <- function() {
  Filter(Negate(is.null), list(
    source    = if (exists("ug_lookup", inherits = TRUE)) ug_lookup else NULL,
    structure = build_cf_lookup(CF_LEVEL),
    moa       = build_moa_lookup()))
}

## ---------------------------------------------------------------------------
## 2. scoring -- all on ONE common node set (present AND labelled by every lens)
## ---------------------------------------------------------------------------
## module-level fit of ONE lens, with a node-label permutation null
lens_modularity <- function(A, lab, nodes, B = LENS_PERM_B, seed = 1) {
  g   <- igraph::graph_from_adjacency_matrix(abs(A[nodes, nodes, drop = FALSE]),
                                             mode = "undirected", weighted = TRUE, diag = FALSE)
  vn  <- igraph::V(g)$name; w <- igraph::E(g)$weight
  grp <- as.integer(factor(lab[vn]))
  Qo  <- igraph::modularity(g, grp, weights = w)
  as  <- tryCatch(suppressWarnings(igraph::assortativity_nominal(g, grp, directed = FALSE)),
                  error = function(e) NA_real_)          # arg name shifted across igraph versions
  set.seed(seed)
  Qn  <- replicate(B, igraph::modularity(g, as.integer(factor(sample(lab[vn]))), weights = w))
  list(Q = Qo, assort = as, k = length(unique(grp)),
       Q_z = if (stats::sd(Qn) > 0) (Qo - mean(Qn)) / stats::sd(Qn) else NA_real_,
       Q_p = (1 + sum(Qn >= Qo)) / (B + 1))
}

## edge-level nested deviance: ONE joint binomial model over all pairs; returns
## each lens's MARGINAL (alone) and UNIQUE (drop-one) deviance + coefficient sign
lens_edge_glm <- function(A, labs, nodes) {
  Aa <- abs(A[nodes, nodes, drop = FALSE]); diag(Aa) <- 0
  ut <- upper.tri(Aa)
  df <- data.frame(y = as.integer(Aa[ut] != 0)); terms <- character(0)
  for (nm in names(labs)) {
    sm <- outer(labs[[nm]][nodes], labs[[nm]][nodes], `==`)
    df[[paste0("same_", nm)]] <- as.integer(sm[ut]); terms <- c(terms, paste0("same_", nm))
  }
  terms <- terms[vapply(terms, function(t) length(unique(df[[t]])) > 1, logical(1))]  # drop degenerate
  if (!length(terms)) return(NULL)
  b0   <- stats::glm(y ~ 1, data = df, family = stats::binomial())
  full <- stats::glm(stats::reformulate(terms, "y"), data = df, family = stats::binomial())
  dev0 <- stats::deviance(b0); devF <- stats::deviance(full)
  do.call(rbind, lapply(terms, function(t) {
    dm <- dev0 - stats::deviance(stats::glm(stats::reformulate(t, "y"), data = df, family = stats::binomial()))
    du <- if (length(terms) > 1)
      stats::deviance(stats::glm(stats::reformulate(setdiff(terms, t), "y"), data = df, family = stats::binomial())) - devF
    else dm
    cf <- unname(coef(full)[t])
    data.frame(lens = sub("^same_", "", t), n_pairs = nrow(df), n_edges = sum(df$y),
               dev_marginal = round(dm, 2), p_marginal = signif(stats::pchisq(dm, 1, lower.tail = FALSE), 3),
               dev_unique = round(du, 2),   p_unique   = signif(stats::pchisq(du, 1, lower.tail = FALSE), 3),
               coef_sign = if (!is.na(cf)) sign(cf) else NA_real_, stringsAsFactors = FALSE)
  }))
}

## score every lens on one unit's network -> tidy per-lens data.frame
lens_compare_unit <- function(A, unit_name, lenses, min_nodes = LENS_MIN_NODES) {
  net <- colnames(A)
  common <- Reduce(intersect, lapply(lenses, function(l) net[!is.na(l[net])]))
  ok <- vapply(lenses, function(l) length(unique(stats::na.omit(l[common]))) >= 2, logical(1))
  lenses <- lenses[ok]
  if (length(common) < min_nodes || length(lenses) < 2) {
    message(sprintf("[lens] %s skipped (common nodes = %d, usable lenses = %d)",
                    unit_name, length(common), length(lenses)))
    return(NULL)
  }
  glm_tab <- lens_edge_glm(A, lenses, common)
  out <- do.call(rbind, lapply(names(lenses), function(nm) {
    m <- lens_modularity(A, lenses[[nm]], common); g <- glm_tab[glm_tab$lens == nm, , drop = FALSE]
    data.frame(unit = unit_name, lens = nm, n_network = length(net), n_common = length(common),
               coverage = round(sum(!is.na(lenses[[nm]][net])) / length(net), 3),
               n_groups = m$k, Q_attr = round(m$Q, 3), Q_z = round(m$Q_z, 2), Q_p = signif(m$Q_p, 3),
               assortativity = round(m$assort, 3),
               dev_marginal = if (nrow(g)) g$dev_marginal else NA_real_,
               dev_unique   = if (nrow(g)) g$dev_unique   else NA_real_,
               p_unique     = if (nrow(g)) g$p_unique     else NA_real_,
               coef_sign    = if (nrow(g)) g$coef_sign    else NA_real_, stringsAsFactors = FALSE)
  }))
  pos <- out[!is.na(out$coef_sign) & out$coef_sign > 0, , drop = FALSE]   # winner: max unique dev, +ve coef
  out$best_lens <- if (nrow(pos)) pos$lens[which.max(pos$dev_unique)] else NA_character_
  out
}

## ---------------------------------------------------------------------------
## 3. driver -- get each unit's glasso network (cache first, else refit), score
## ---------------------------------------------------------------------------
## Uses load_fitted() when a cached fitted_*.rds exists; otherwise refits with
## prep_unit()+fit_glasso() (same functions, same defaults as the mainline).
lens_get_A <- function(u, kind, rep.num = STARS_REPNUM) {
  fit <- load_fitted(kind, u$safe)
  if (!is.null(fit) && !is.null(fit$A_glasso)) return(fit$A_glasso)
  pu <- tryCatch(prep_unit(kind, u$site_ids), error = function(e) NULL)
  if (is.null(pu) || pu$n_comp < 5 || pu$n_sites < 4) return(NULL)
  fit_glasso(pu$Xlog, rep.num = rep.num)$A
}

lens_run <- function(unit_list, kind = "conc", lenses = NULL, rep.num = STARS_REPNUM,
                     out_dir = file.path(outroot, "lenses")) {
  if (is.null(lenses)) lenses <- assemble_lenses()
  message("[lens] lenses: ", paste(names(lenses), collapse = ", "))
  if (length(lenses) < 2) { message("[lens] need >= 2 lenses (run ClassyFire chunk / add compound_moa.csv)"); return(NULL) }
  master <- do.call(rbind, Filter(Negate(is.null), lapply(unit_list, function(u) {
    A <- tryCatch(lens_get_A(u, kind, rep.num), error = function(e) NULL)
    if (is.null(A)) return(NULL)
    tryCatch(lens_compare_unit(A, u$safe, lenses), error = function(e) { message("[lens] ", u$safe, ": ", conditionMessage(e)); NULL })
  })))
  if (is.null(master)) return(NULL)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(master, file.path(out_dir, "master_lens_comparison.csv"), row.names = FALSE)
  verdict <- do.call(rbind, by(master, master$unit, function(d) data.frame(
    unit = d$unit[1], best_lens = d$best_lens[1],
    dev_unique_source    = d$dev_unique[d$lens == "source"][1],
    dev_unique_structure = d$dev_unique[d$lens == "structure"][1],
    dev_unique_moa       = if ("moa" %in% d$lens) d$dev_unique[d$lens == "moa"][1] else NA_real_,
    n_common = d$n_common[1], stringsAsFactors = FALSE)))
  utils::write.csv(verdict, file.path(out_dir, "master_lens_verdict.csv"), row.names = FALSE)
  list(master = master, verdict = verdict, out_dir = out_dir)
}

## two-panel summary figure: unique deviance (A) + modularity-z (B), by lens/unit
lens_figure <- function(master, out_dir = file.path(outroot, "lenses")) {
  cols <- c(source = "#457B9D", structure = "#2A9D8F", moa = "#E9C46A")
  units_ord <- unique(master$unit); lenses <- intersect(names(cols), unique(master$lens))
  as_mat <- function(col) { m <- matrix(NA_real_, length(lenses), length(units_ord), dimnames = list(lenses, units_ord))
    for (i in seq_len(nrow(master))) m[master$lens[i], master$unit[i]] <- master[[col]][i]; m }
  M_dev <- as_mat("dev_unique"); M_z <- as_mat("Q_z"); M_sg <- as_mat("coef_sign")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  grDevices::png(file.path(out_dir, "fig_lens_comparison.png"), width = 1600, height = 720, res = 150)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mfrow = c(1, 2), mar = c(9, 4.5, 3, 1))
  bp <- graphics::barplot(M_dev, beside = TRUE, col = cols[lenses], border = NA, las = 2,
                          ylab = "Unique deviance (edge model, 1 df)", main = "A  What each lens adds over the others")
  neg <- which(!is.na(M_sg) & M_sg < 0, arr.ind = TRUE)                 # hatch negative-coefficient bars
  if (length(neg)) for (r in seq_len(nrow(neg)))
    graphics::rect(bp[neg[r,1], neg[r,2]] - .5, 0, bp[neg[r,1], neg[r,2]] + .5, M_dev[neg[r,1], neg[r,2]],
                   density = 18, col = "grey30", border = NA)
  graphics::abline(h = stats::qchisq(0.95, 1), lty = 3, col = "grey40")
  graphics::legend("topright", legend = lenses, fill = cols[lenses], border = NA, bty = "n")
  graphics::barplot(M_z, beside = TRUE, col = cols[lenses], border = NA, las = 2,
                    ylab = "Modularity z (label-permutation null)", main = "B  Do same-label compounds form the modules?")
  graphics::abline(h = 1.96, lty = 3, col = "grey40")
  invisible(file.path(out_dir, "fig_lens_comparison.png"))
}
