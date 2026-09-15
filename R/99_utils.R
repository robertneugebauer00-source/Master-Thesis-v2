## =============================================================================
## 99_utils.R -- figure refresh and HTML report builder (no analysis)
## Split out of Clean_Workflow_vRN-14.08.26.Rmd on 14.09.2026 -- code is verbatim.
## Original chunks: vrn-refresh-figures | build-report
## [repo 15.09.2026] no longer verbatim: sourcing defines functions only
## (refresh_all_figures / build_review_report / main); refresh reuses cached
## fits; reach_units moved to 01_data.R; metric_row() moved to 02_functions.R.
## =============================================================================

## ---- vrn-refresh-figures  [original L2097-2174] ----
## ---------------------------------------------------------------------------
## [vRN 14.08.2026] figures-only refresh: fig1 + fig1b for every unit.
## Skips null_test() and the WGCNA diagnostic/dashboard figures on purpose.
##
## [repo 15.09.2026] NOW ACTUALLY "WITHOUT REFITTING". The mainline run saves
## fitted_<kind>_<unit>.rds per unit (adjacency matrices + memberships); when
## that artefact exists we redraw from it and skip fit_glasso/fit_mb/fit_wgcna
## entirely. Only units with no cached fit fall back to refitting (and get
## refit = TRUE in the log). Note the cached fit reflects the parameters of the
## run that produced it -- if you changed PREV_THRESH etc., re-run the mainline
## instead of refreshing.
## ---------------------------------------------------------------------------

refresh_unit_figures <- function(kind, unit, rep.num = STARS_REPNUM) {
  kind_label <- if (kind == "conc") "Concentrations (ng/L)" else "Toxic units (TRIDENT EC10)"
  outdir <- file.path(outroot, kind_dir(kind), unit$safe)
  if (!dir.exists(outdir)) return(NULL)          # unit was never run -> nothing to refresh

  fc <- load_fitted(kind, unit$safe)
  ## a glasso-only artefact (e.g. written by the optional _targets.R pipeline,
  ## which runs the structural core but not MB/WGCNA) cannot feed the
  ## three-method figure panels -- treat it as no cache and refit
  if (!is.null(fc) && (is.null(fc$A_mb) || is.null(fc$A_wgcna))) fc <- NULL
  if (!is.null(fc)) {
    gl <- list(A = fc$A_glasso); mb <- list(A = fc$A_mb)
    wg <- list(A = fc$A_wgcna, modules = fc$mods_wgcna)
    mods_gl <- fc$mods_glasso; mods_mb <- fc$mods_mb
    n_sites <- fc$n_sites; n_comp <- fc$n_comp
  } else {
    pu <- prep_unit(kind, unit$site_ids)
    if (pu$n_comp < 5 || pu$n_sites < 4) return(NULL)
    gl <- fit_glasso(pu$Xlog, rep.num = rep.num)
    mb <- fit_mb(pu$Xlog,     rep.num = rep.num)
    wg <- fit_wgcna(pu$Xlog)
    mods_gl <- greedy_modules(gl$A); mods_mb <- greedy_modules(mb$A)
    n_sites <- pu$n_sites; n_comp <- pu$n_comp
  }

  ## [repo 15.09.2026] metric_row() from 02_functions.R -- the local copy that
  ## used to live here had drifted from the one in run_unit() (no Nodes (pool)).
  metrics <- list(Glasso = metric_row(metrics_from_adj(gl$A, mods_gl)),
                  MB     = metric_row(metrics_from_adj(mb$A, mods_mb)),
                  WGCNA  = metric_row(metrics_from_adj(wg$A, wg$modules)))

  rd <- tryCatch(risk_driver_table(unit$site_ids, gl$A), error = function(e) NULL)
  rd_names <- if (is.null(rd)) character(0) else rd$Compound[rd$absolute_risk_driver %in% TRUE]

  As     <- list(Glasso = gl$A, MB = mb$A, WGCNA = wg$A)
  by_mod <- list(Glasso = node_cols_module(gl$A, mods_gl),
                 MB     = node_cols_module(mb$A, mods_mb),
                 WGCNA  = node_cols_module(wg$A, wg$modules))
  by_ug  <- list(Glasso = node_cols_ug(gl$A), MB = node_cols_ug(mb$A), WGCNA = node_cols_ug(wg$A))
  mods   <- list(Glasso = mods_gl, MB = mods_mb, WGCNA = wg$modules)

  plot_modules(As, by_mod, mods, metrics, file.path(outdir, "fig1_network_panel.png"),
               unit$label, kind_label, risk = rd_names)
  plot_usegroups(As, by_ug, metrics, file.path(outdir, "fig1b_network_panel_usegroup.png"),
                 unit$label, kind_label, ug_present_in(As), risk = rd_names)
  ## fig4 carries the same driver marker, so it is stale too; rd is already in hand
  if (!is.null(rd)) plot_hub_vs_risk(rd, file.path(outdir, "fig4_hub_vs_riskdriver.png"), unit$label)

  data.frame(kind = kind, unit = unit$safe, n_sites = n_sites, n_comp = n_comp,
             n_drivers = length(rd_names), refit = is.null(fc), row.names = NULL)
}

## [repo 15.09.2026] reach_units now comes from 01_data.R (single source of
## truth; the local copy that used to live here was removed).

## [repo 15.09.2026] the refresh loop is a FUNCTION now. Sourcing this file
## defines functions and builds NOTHING; call refresh_all_figures() explicitly
## (run_all.R does this for the "utils" step).
refresh_all_figures <- function(kinds = c("conc", "tu"), rep.num = STARS_REPNUM) {
  refresh_log <- list()
  for (kind in kinds) {
    for (u in c(units, basin_units, section_units, reach_units)) {
      r <- tryCatch(refresh_unit_figures(kind, u, rep.num = rep.num),
                    error = function(e) { message("[vRN] refresh failed: ", kind, " ", u$safe,
                                                  " -- ", conditionMessage(e)); NULL })
      if (!is.null(r)) { refresh_log[[length(refresh_log) + 1]] <- r; message("[vRN] refreshed ", kind, " ", u$safe) }
    }
  }
  refresh_log <- if (length(refresh_log)) do.call(rbind, refresh_log) else NULL
  print(refresh_log)
  invisible(refresh_log)
}

## ---- build-report  [original L2593-2713] ----

build_review_report <- function(outroot, out_html = file.path(outroot, "REVIEW_report.html")) {

  esc <- function(x) { x <- as.character(x)
    x <- gsub("&","&amp;",x,fixed=TRUE); x <- gsub("<","&lt;",x,fixed=TRUE)
    gsub(">","&gt;",x,fixed=TRUE) }

  img_tag <- function(path) {
    if (!file.exists(path)) return("<p style='color:#999'>[figure missing]</p>")
    uri <- base64enc::dataURI(file = path, mime = "image/png")
    sprintf('<img loading="lazy" src="%s" style="max-width:100%%;height:auto;border:1px solid #e0e0e0;border-radius:6px;margin:8px 0;">', uri)
  }

  tbl_tag <- function(path, id) {
    if (!file.exists(path)) return("")
    df <- tryCatch(read.csv(path, check.names = FALSE), error = function(e) NULL)
    if (is.null(df) || !nrow(df)) return("")
    th <- paste0("<th>", esc(colnames(df)), "</th>", collapse = "")
    body <- apply(df, 1, function(r) paste0("<tr>", paste0("<td>", esc(r), "</td>", collapse = ""), "</tr>"))
    sprintf('<table id="%s" class="review display compact" style="width:100%%"><thead><tr>%s</tr></thead><tbody>%s</tbody></table>',
            id, th, paste(body, collapse = ""))
  }

  section_unit <- function(kind_dir, safe, label) {
    d <- file.path(outroot, kind_dir, safe)
    if (!dir.exists(d)) return("")
    idp <- gsub("[^A-Za-z0-9]+", "_", paste(kind_dir, safe, sep = "_"))
    paste0(
      sprintf('<h3 id="%s">%s</h3>', idp, esc(label)),
      '<div class="grid2">',
        '<div>', '<h4>Networks — modules</h4>', img_tag(file.path(d, "fig1_network_panel.png")), '</div>',
        '<div>', '<h4>Networks — use groups</h4>', img_tag(file.path(d, "fig1b_network_panel_usegroup.png")), '</div>',
      '</div>',
      '<h4>Sentinel vs. absolute risk driver</h4>', img_tag(file.path(d, "fig4_hub_vs_riskdriver.png")),
      '<div class="grid2">',
        '<div>', '<h4>WGCNA diagnostics</h4>', img_tag(file.path(d, "fig2_wgcna_diagnostics.png")), '</div>',
        '<div>', '<h4>WGCNA dashboard</h4>', img_tag(file.path(d, "fig3_wgcna_dashboard.png")), '</div>',
      '</div>',
      '<h4>Topology</h4>',            tbl_tag(file.path(d, "table1_topology_summary.csv"),   paste0("t1_", idp)),
      '<h4>Module ARI</h4>',          tbl_tag(file.path(d, "table2b_ARI_matrix.csv"),        paste0("t2_", idp)),
      '<h4>Null-model test</h4>',    tbl_tag(file.path(d, "table3_nullmodel_sensitivity.csv"), paste0("t3_", idp)),
      '<h4>Absolute risk drivers</h4>', tbl_tag(file.path(d, "table4_risk_drivers.csv"),    paste0("t4_", idp)),
      '<hr>'
    )
  }

  # collect units present on disk, ordered group-then-basin
  present <- list.dirs(file.path(outroot, "concentration"), recursive = FALSE, full.names = FALSE)
  lab_of  <- function(safe) gsub("_", " ", safe)

  build_kind <- function(kind_dir, title) {
    us <- list.dirs(file.path(outroot, kind_dir), recursive = FALSE, full.names = FALSE)
    if (!length(us)) return("")
    paste0(sprintf('<h2 id="sec_%s">%s</h2>', kind_dir, esc(title)),
           paste0(vapply(us, function(s) section_unit(kind_dir, s, lab_of(s)), character(1)), collapse = ""))
  }

  master <- function() {
    m <- c("master_topology_summary.csv","master_ARI.csv","master_nullmodel_sensitivity.csv","master_risk_drivers.csv")
    lab <- c("Topology (all units)","Module ARI (all units)","Null-model (all units)","Absolute risk drivers (all units)")
    paste0('<h2 id="sec_master">Master summary tables</h2>',
           paste0(mapply(function(f, l) paste0(sprintf("<h4>%s</h4>", l),
                    tbl_tag(file.path(outroot, f), paste0("m_", gsub("[^A-Za-z0-9]+","_", f)))),
                    m, lab), collapse = ""))
  }

  ## [vPI2] cross-basin topology comparison — flagged as a known gap
  basin_topo_section <- function() {
    f <- file.path(outroot, "master_topology_summary.csv")
    if (!file.exists(f)) return("")
    d <- read.csv(f, check.names = FALSE)
    d <- d[d$Data == "Concentration" & d$Scale == "basin" & d$Model == "Glasso",
           c("Unit","Sites","Nodes","Edges","Density","Modularity_Q","Avg_path_length")]
    if (!nrow(d)) return("")
    d <- d[order(-d$Sites), ]
    th   <- paste0("<th>", esc(colnames(d)), "</th>", collapse = "")
    body <- apply(d, 1, function(r) paste0("<tr>", paste0("<td>", esc(r), "</td>", collapse = ""), "</tr>"))
    tab  <- sprintf('<table class="review display compact" style="width:100%%"><thead><tr>%s</tr></thead><tbody>%s</tbody></table>',
                    th, paste(body, collapse = ""))
    paste0(
      '<h2 id="sec_basin_topo">Cross-basin topology — a known gap</h2>',
      '<div class="note"><b>This is a gap, not a finished result.</b> The per-basin topology metrics below are computed and available, but there is currently <b>no explicit cross-basin comparison</b> (no test, no module-sharing analysis), and the raw comparison is <b>confounded by sample size</b>: n ranges from ~20 to ~150 sites. Sparser networks are mechanically more fragmented and more modular, so a basin with few sites (e.g. Basque, Tajo) will look more modular for reasons unrelated to its chemistry. To compare basins fairly, all basins should be <b>rarefied to a common site count</b> before their topologies are compared, and a cross-basin module-sharing check added (do the same compounds cluster together across basins?). Treat the table below as raw input to that future analysis.</div>',
      tab)
  }

  html <- paste0(
    '<!DOCTYPE html><html><head><meta charset="utf-8"><title>Network review — vRN2</title>',
    '<link rel="stylesheet" href="https://cdn.datatables.net/1.13.6/css/jquery.dataTables.min.css">',
    '<style>',
      'body{font-family:-apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif;margin:0;color:#222;line-height:1.5}',
      '.wrap{max-width:1200px;margin:0 auto;padding:24px}',
      'h1{color:#264653} h2{color:#2A9D8F;border-bottom:2px solid #e9ecef;padding-bottom:6px;margin-top:40px}',
      'h3{color:#264653;margin-top:28px} h4{color:#555;margin:14px 0 4px}',
      '.grid2{display:grid;grid-template-columns:1fr 1fr;gap:16px}',
      '@media(max-width:900px){.grid2{grid-template-columns:1fr}}',
      'table.review{font-size:12px;border-collapse:collapse} table.review td,table.review th{border:1px solid #e6e6e6;padding:3px 6px}',
      '.toc{background:#f8f9fa;border:1px solid #e9ecef;border-radius:8px;padding:12px 18px;margin:18px 0}',
      '.note{background:#fff8e1;border-left:4px solid #E9C46A;padding:10px 14px;border-radius:4px;margin:14px 0;font-size:14px}',
    '</style></head><body><div class="wrap">',
    '<h1>Chemical co-occurrence networks — vPI3 review</h1>',
    sprintf('<p style="color:#777">Generated %s · glasso (primary) · MB (cross-check) · WGCNA (modules) · absolute-risk-driver overlay (RQ ≥ %.2f).</p>', Sys.Date(), RQ_THRESH),
    '<div class="note"><b>How to read this:</b> violet rings (network panels) and red triangles (scatter) mark <b>absolute risk drivers</b> (RQ ≥ 0.02 at ≥1 site). Black rings mark <b>structural hubs (sentinels)</b>. A risk driver that sits far left in the sentinel-vs-driver plot (low centrality) is a <b>peripheral / off-network</b> risk the co-occurrence network alone would miss.</div>',
    '<div class="toc"><b>Contents</b><br><a href="#sec_master">Master tables</a> · <a href="#sec_basin_topo">Cross-basin topology</a> · <a href="#sec_concentration">Concentration networks</a> · <a href="#sec_toxic_unit">Toxic-unit networks</a></div>',
    master(),
    basin_topo_section(),
    build_kind("concentration", "Concentration networks (structure &amp; sources)"),
    build_kind("toxic_unit",    "Toxic-unit networks (near-identical wiring; kept for reference)"),
    '<script src="https://code.jquery.com/jquery-3.7.1.min.js"></script>',
    '<script src="https://cdn.datatables.net/1.13.6/js/jquery.dataTables.min.js"></script>',
    '<script>$(function(){$("table.review").each(function(){$(this).DataTable({paging:$(this).find("tbody tr").length>15,searching:true,info:false,order:[]});});});</script>',
    '</div></body></html>'
  )
  writeLines(html, out_html, useBytes = TRUE)
  message("Report written: ", out_html)
  invisible(out_html)
}

## [repo 15.09.2026] entry point. Sourcing this file no longer runs anything;
## run_all.R calls main() for its "utils" step, and `Rscript R/99_utils.R`
## triggers it via the sys.nframe() guard below.
main <- function() {
  refresh_all_figures()
  build_review_report(outroot)
}
if (sys.nframe() == 0L) main()


