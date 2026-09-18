## =============================================================================
## 01_data.R -- load data, MDL floor, toxic units, spatial unit definitions
## Split out of Clean_Workflow_vRN-14.08.26.Rmd on 14.09.2026 -- code is verbatim.
## Original chunks: #DATA LOADING v1.0 | mdl-lookup | #TOXIC UNIT CALC | groups
## =============================================================================

## ---- #DATA LOADING v1.0  [original L427-463] ----

# --- site metadata (sheet B) ---

b_raw <- read_excel(conc_file, sheet = "B_parameters_info", col_names = FALSE)
hdr   <- which(b_raw[[1]] == "Sample code")

sites_meta <- b_raw[(hdr + 1):nrow(b_raw), ]
colnames(sites_meta) <- as.character(b_raw[hdr, ])
sites_meta <- sites_meta %>%
  rename(site_id     = `Sample code`,
         lat         = Latitude,
         lon         = Longitude,
         country     = Country,
         river_basin = `River basin`) %>%
  mutate(lat = as.numeric(lat), lon = as.numeric(lon))

# --- concentration matrix (compounds in rows, sites in columns) ---

ab       <- read_excel(conc_file, sheet = "A_B_results", col_names = FALSE)
name_row <- which(ab[[1]] == "Name")
first_eus <- which(grepl("^EUS_", as.character(ab[name_row, ])))[1]
data_rows <- (name_row + 1):nrow(ab)

compounds  <- as.character(ab[data_rows, 1, drop = TRUE])
use_groups <- as.character(ab[data_rows, 2, drop = TRUE])
compound_inchikey <- as.character(ab[data_rows, 6, drop = TRUE])   ## [vPI3] "InChIKey" column, used to match MDLs

conc_mat <- apply(ab[data_rows, first_eus:ncol(ab)], 2, as.numeric)
rownames(conc_mat) <- compounds
colnames(conc_mat) <- as.character(ab[name_row, first_eus:ncol(ab)])
conc_mat[is.na(conc_mat)] <- 0          # non-detect = 0 (below MDL)

# compound -> use group, for colouring later

ug_lookup <- setNames(use_groups, compounds)

## ---- mdl-lookup  [original L547-580] ----

## [vRN] guarded: if the MDL file is not on disk yet (still on Marvin), fall back to the
## +DEFAULT_PSEUDO_CONC ng/L floor for every compound instead of erroring on read.csv().
if (file.exists(mdl_file)) {
  mdl_tab <- read.csv(mdl_file, stringsAsFactors = FALSE)

  mdl_by_ik   <- setNames(mdl_tab$mdl_value_ng_L, toupper(trimws(mdl_tab$inchikey)))
  mdl_by_name <- setNames(mdl_tab$mdl_value_ng_L, tolower(trimws(mdl_tab$chemical_name)))

  ik_key   <- toupper(trimws(compound_inchikey))
  name_key <- tolower(trimws(compounds))

  mdl_ngL <- setNames(mdl_by_ik[ik_key], compounds)
  unmatched <- is.na(mdl_ngL)
  mdl_ngL[unmatched] <- mdl_by_name[name_key[unmatched]]
  unmatched <- is.na(mdl_ngL)

  message("MDL match: ", sum(!unmatched), "/", length(compounds),
          " compounds matched to MDL_CECs_summary.csv; ", sum(unmatched),
          " fall back to the default +", DEFAULT_PSEUDO_CONC, " ng/L pseudocount: ",
          paste(compounds[unmatched], collapse = ", "))

  mdl_half_ngL <- mdl_ngL / 2
  mdl_half_ngL[unmatched] <- DEFAULT_PSEUDO_CONC   # same fallback value as vPI2, now per compound
} else {                                                    ## [vRN] MDL file absent
  message("[vRN] MDL file not found at ", mdl_file,
          " -- download MDL_CECs_summary.csv from Marvin into raw_data/. ",
          "Falling back to +", DEFAULT_PSEUDO_CONC, " ng/L for ALL compounds (vPI2 behaviour).")
  mdl_ngL      <- setNames(rep(NA_real_, length(compounds)), compounds)
  mdl_half_ngL <- setNames(rep(DEFAULT_PSEUDO_CONC, length(compounds)), compounds)
  unmatched    <- rep(TRUE, length(compounds))
}

## ---- #TOXIC UNIT CALC  [original L585-619] ----

## [vRN] robust EC10 loader: raw TRIDENT (pred + map) if present, else the processed
## TRIDENT_EC10_ngL.csv already written to data_for_marvin, else skip TU gracefully.
if (file.exists(pred_csv) && file.exists(map_csv)) {
  pred <- read.csv(pred_csv, check.names = FALSE)
  map  <- read.csv(map_csv,  check.names = FALSE)
  ec10_ngL <- map %>%
    select(cmpdname, SMILES) %>%
    left_join(pred %>% transmute(SMILES, mg = `predictions (mg/L)`), by = "SMILES") %>%
    { setNames(.$mg * 1e6, .$cmpdname) }
} else if (file.exists(ec10_fallback)) {
  message("[vRN] raw TRIDENT files not found -> using processed EC10 from ", basename(ec10_fallback))
  ef <- read.csv(ec10_fallback, check.names = FALSE)
  ec10_ngL <- setNames(as.numeric(ef$EC10_ngL), ef$Compound)
} else {
  warning("[vRN] no EC10 source found -> toxic-unit (TU) analysis will be skipped.")
  ec10_ngL <- numeric(0)
}
ec10_ngL <- ec10_ngL[is.finite(ec10_ngL)]

# Compounds that have an EC10 can get a toxic unit  (TU = MEC / EC10 = per-site RQ)

have_ec10 <- intersect(rownames(conc_mat), names(ec10_ngL))
## @KEY  toxic unit TU = MEC / EC10. Risk scale, NOT a second network ---------------------------------
tu_mat    <- sweep(conc_mat[have_ec10, ], 1, ec10_ngL[have_ec10], "/")
tu_mat[!is.finite(tu_mat)] <- 0

## [vPI3] same MDL/2 floor, converted to toxic-unit scale (TU = MEC / EC10) for
## compounds that have both an MDL and a TRIDENT EC10. Compounds missing either
## fall back to prep_unit()'s vPI2-style within-unit 0.5*min(nonzero) floor.
mdl_half_tu <- (mdl_half_ngL[have_ec10] / ec10_ngL[have_ec10])
mdl_half_tu <- mdl_half_tu[is.finite(mdl_half_tu)]


## ---- groups  [original L626-723] ----
groups <- list(
  list(label = "Upper Danube",  safe = "Danube_Upper",  basin = "Danube",
       countries = c("Germany","Austria","Czechia")),
  list(label = "Middle Danube", safe = "Danube_Middle", basin = "Danube",
       countries = c("Slovenia","Slovakia","Hungary","Croatia")),
  list(label = "Lower Danube",  safe = "Danube_Lower",  basin = "Danube",
       countries = c("Serbia","Bulgaria","Romania","Moldova","Ukraine")),
  list(label = "Upper Rhine",   safe = "Rhine_Upper",   basin = "Rhine",
       countries = c("Switzerland")),
  list(label = "Lower Rhine",   safe = "Rhine_Lower",   basin = "Rhine",
       countries = c("Germany","Netherlands"))
)

# Attach the site IDs for each group

groups <- lapply(groups, function(g) {
  g$site_ids <- sites_meta %>%
    filter(river_basin == g$basin, country %in% g$countries) %>%
    pull(site_id)
  g
})

# Note: no pooled "Rhine+Danube" unit here. Each basin is run independently
# as its own whole-basin unit further below (see `basin_units`).
## [vRN 03.09.2026] Apply the same minimum-n floor to the country groups that already
## applies to basins and sections. Without it, "Upper Rhine" = Switzerland = 13 sites was
## carried as a unit and returned ZERO edges in every run -- and cannot be rescued: 0 edges
## at MIN_DETECT_UNIT 2, 3, 4, 5 and 6 alike (outputs_vRN/rhine_guard_sweep_20260903.csv).
## The Rhine is still covered, by Basin_Rhine and by the two river sections below.
units <- Filter(function(g) length(g$site_ids) >= BASIN_MIN_N, groups)
message("[vRN] country groups (n >= ", BASIN_MIN_N, "): ",
        paste(sprintf("%s(%d)", sapply(units, `[[`, "safe"),
                      sapply(units, function(u) length(u$site_ids))), collapse = ", "),
        if (length(units) < length(groups))
          paste0("  [dropped: ",
                 paste(setdiff(sapply(groups, `[[`, "safe"), sapply(units, `[[`, "safe")),
                       collapse = ", "), "]") else "")

## [vPI2] ---- river-basin-specific units: all sites of a basin, n >= BASIN_MIN_N ----
basin_counts <- sort(table(sites_meta$river_basin), decreasing = TRUE)
big_basins   <- names(basin_counts)[basin_counts >= BASIN_MIN_N]
basin_units  <- lapply(big_basins, function(bn) {
  list(label   = paste0(bn, " (whole basin)"),
       safe    = paste0("Basin_", gsub("[^A-Za-z0-9]+", "_", bn)),
       basin   = bn,
       site_ids = sites_meta$site_id[sites_meta$river_basin == bn])
})
message("Basin-specific units (n >= ", BASIN_MIN_N, "): ",
        paste(sprintf("%s(%d)", big_basins, basin_counts[big_basins]), collapse = ", "))

## [vRN] ---- latitude/longitude Upper/Middle/Lower sections (incl. the Elbe) ----
## Elbe & Rhine flow roughly S->N, so latitude tracks the river axis (upstream = low lat);
## the Danube flows W->E, so longitude is the right axis (upstream = low lon). We split each
## basin's sites into along-river thirds and label them so "Upper" is always upstream.
SECTION_SPEC <- list(
  Elbe   = list(axis = "lat", upstream = "low"),
  Rhine  = list(axis = "lat", upstream = "low"),
  Danube = list(axis = "lon", upstream = "low")
)
## [vRN 03.09.2026] The number of sections is DERIVED from basin size, not fixed at 3.
## Evidence (outputs_vRN/rhine_resplit_20260903.csv + rhine_guard_sweep_20260903.csv):
##   Rhine thirds = 21 sites each. The middle third returns ZERO edges -- StARS picks the
##   largest lambda on the grid because the first non-empty graph already sits at
##   instability 0.063 against a 0.05 threshold. Raising MIN_DETECT_UNIT does not rescue
##   it (z = 0.22 at guard 4, -0.80 at guard 5) and costs signal everywhere else.
##   Rhine halves = 30 sites each: both return dense, significant networks.
## So a section needs SECTION_MIN_N sites, and a basin gets floor(n / SECTION_MIN_N)
## sections, capped at 3. Rhine -> 2, Elbe and Danube -> 3.
make_sections <- function(bn, spec) {
  ax <- spec$axis
  sm <- sites_meta[sites_meta$river_basin == bn & is.finite(sites_meta[[ax]]), ]
  k  <- min(3L, as.integer(floor(nrow(sm) / SECTION_MIN_N)))
  if (k < 2L) return(list())                    # too few sites to split this basin at all
  v    <- sm[[ax]]
  qs   <- quantile(v, seq_len(k - 1L) / k, na.rm = TRUE)
  tile <- cut(v, breaks = c(-Inf, qs, Inf), labels = paste0("t", seq_len(k)))
  pos_names <- if (k == 2L) c("Upper", "Lower") else c("Upper", "Middle", "Lower")
  posmap <- setNames(pos_names, paste0("t", seq_len(k)))
  if (spec$upstream != "low") posmap <- setNames(rev(pos_names), paste0("t", seq_len(k)))
  lapply(pos_names, function(pos) {
    key <- names(posmap)[posmap == pos]
    list(label    = paste(pos, bn, "(section)"),
         safe     = paste0("Sec_", gsub("[^A-Za-z0-9]+", "_", bn), "_", pos),
         basin    = bn,
         site_ids = sm$site_id[tile == key])
  })
}
section_units <- unlist(lapply(intersect(names(SECTION_SPEC), big_basins),
                               function(bn) make_sections(bn, SECTION_SPEC[[bn]])),
                        recursive = FALSE)
section_units <- Filter(function(u) length(u$site_ids) >= BASIN_MIN_N, section_units)
message("[vRN] river sections (n >= ", BASIN_MIN_N, "): ",
        if (length(section_units))
          paste(sprintf("%s(%d)", sapply(section_units, `[[`, "label"),
                        sapply(section_units, function(u) length(u$site_ids))), collapse = ", ")
        else "none")

## [repo 15.09.2026] ---- upstream/downstream REACH units ----------------------
## SINGLE SOURCE OF TRUTH for what a reach is. This definition used to live in
## vRN_followups.R (analysis A4a, outside this repo) with a second copy in
## 99_utils.R -- two sources of truth that could drift. It now lives here, next
## to the other unit definitions; 99_utils.R and vRN_followups.R should both
## use `reach_units` from this file.
make_reach <- function(basin, axis) {
  sm <- as.data.frame(sites_meta)
  d  <- sm[!is.na(sm$river_basin) & sm$river_basin == basin & is.finite(sm[[axis]]), ]
  d  <- d[d$site_id %in% colnames(conc_mat), ]
  m  <- median(d[[axis]])
  list(list(label = paste(basin, "- Upstream"),   safe = paste0("Reach_", basin, "_Upstream"),
            basin = basin, site_ids = d$site_id[d[[axis]] <  m]),
       list(label = paste(basin, "- Downstream"), safe = paste0("Reach_", basin, "_Downstream"),
            basin = basin, site_ids = d$site_id[d[[axis]] >= m]))
}
## same median split as the original vRN_followups.R definition (A4a)
reach_units <- c(make_reach("Danube", "lon"), make_reach("Elbe", "lat"))

