## =============================================================================
## 00_setup.R -- palette, packages, file paths, ALL analysis parameters
## Split out of Clean_Workflow_vRN-14.08.26.Rmd on 14.09.2026 -- code is verbatim.
## Original chunks: #SETUP FOR PALETTE etc. | #PRACKAGE AND FILE SETUP
## =============================================================================

## ---- #SETUP FOR PALETTE etc.  [original L260-289] ----
# ===========================================================================
# Shared setup - compound use-group colour palette
# Setup for every Graph ( for Consistency)
# ===========================================================================

## rm(list = ls())   ## [split 14.09.2026] DISABLED -- harmless in a knit, destructive when sourced into a live session

ug_palette <- c(
  "pharmaceutical"      = "#E63946",
  "pesticide"           = "#2A9D8F",
  "sweetener"           = "#E9C46A",
  "industrial"          = "#457B9D",
  "corrosion inhibitor" = "#9B2226",
  "human metabolite"    = "#F4A261",
  "UV filter"           = "#A8DADC",
  "biocide"             = "#264653",
  "stimulans"           = "#BFAE48",
  "repellent"           = "#6A994E",
  "flame retardant"     = "#BC6C25",
  "surfactant"          = "#7678ED",
  "food ingredient"     = "#F9C784",
  "dye"                 = "#C77DFF",
  "PFC"                 = "#48CAE4",
  "plastic additive"    = "#B5838D",
  "rubber additive"     = "#6B705C",
  "bittern"             = "#A5A58D",
  "other"               = "#BBBBBB"
)

## ---- #PRACKAGE AND FILE SETUP  [original L293-421] ----

# install.packages(c("readxl","tidyverse","igraph","huge","magick","openxlsx",
#                    "base64enc"))                    ## [vPI2] base64enc for the report
# BiocManager::install("WGCNA")

library(readxl)
library(tidyverse)
library(igraph)
library(huge)
library(WGCNA)
library(dynamicTreeCut)
library(magick)
library(uwot)
library(openxlsx)
library(base64enc)      ## [vPI2] embed figures in the standalone HTML report
if (!requireNamespace("pagedown", quietly = TRUE)) message("[vRN] pagedown not installed — skipping (only used for optional PDF export of the report).") else library(pagedown)

#File Setup

## [repo 14.09.2026] Data root. Set SGH_PROJ to point at your own copy of the data
## tree (raw_data/, outputs/, scripts/); the default is the original Windows path.
proj      <- Sys.getenv("SGH_PROJ", unset = "C:/Users/rober/OneDrive/Masterkram/R")
conc_file <- file.path(proj, "scripts", "Finckh_Carmona_2023_PANGAEA_R1.xlsx")     ## [vRN]

# TRIDENT EC10 predictions + the compound -> SMILES map used to join them

pred_csv  <- file.path(proj, "raw_data", "TRIDENT_prediction_results_EC10_GRO.csv")   ## [vRN] raw TRIDENT (optional)
map_csv       <- file.path(proj, "outputs", "danube conc vs TU", "trident", "Danube_TRIDENT_mapping.csv")  ## [vRN]
ec10_fallback <- file.path(proj, "outputs_updated", "data_for_marvin", "TRIDENT_EC10_ngL.csv")             ## [vRN] processed EC10 fallback

outroot   <- file.path(proj, "outputs_vRN")   ## [vRN] separate output tree (Pedro's vPI3 tree left intact)

## [vPI3] method detection limits, same lab/method across projects, different years
## (see raw_data/MDL_CECs_summary.csv, one row per compound, lowest reported MDL across sources)
mdl_file <- file.path(proj, "raw_data", "MDL_CECs_summary.csv")   ## [vRN] download Pedro's Marvin file here
DEFAULT_PSEUDO_CONC <- 1   # ng/L fallback for compounds with no MDL match (kept from vPI2 for continuity)

## [vPI2] global analysis parameters (were hard-coded before)
## [vRN 05.09.2026] PREVALENCE FILTER -- PROVENANCE AND PURPOSE (meeting TODO 04.09.2026).
##   ORIGINAL REFERENCE. The 10 % comes from Hernandez et al. (2021, ISME J 15:1722),
##     the SGH microbial-network paper: "co-occurrence networks using SparCC ... with
##     OTUs present in >= 10 % of all samples [17]". Their [17] is Herren & McMahon
##     (2017, ISME J 11:2426), which is where the practice actually originates and which
##     uses 5 %: "we removed taxa that were not present in at least 5 % of samples, as we
##     were not confident that we could recover robust connectedness estimates for very
##     rare taxa". Data type there: phytoplankton counts plus one 16S amplicon set.
##     So the stated purpose is ESTIMATOR RELIABILITY for rare variables -- not
##     compositionality, not sparsification, and the 10 % is a doubled convention with
##     no derivation behind it.
##   IS IT OPPOSED TO THE GLASSO SPARSIFICATION? No -- the two act on different objects.
##     The filter selects the NODE SET (which variables enter the model at all); glasso +
##     StARS select the EDGE SET (which conditional dependencies are non-zero given those
##     nodes). It is not double sparsification of the same quantity.
##   WHAT TRANSFERS TO CHEMICAL DATA. The microbiome rationale (rounded zeros, sequencing
##     depth, compositional closure) does NOT transfer: here a zero is a true
##     non-detection, i.e. LEFT CENSORING. What does transfer is the support argument --
##     a mostly-censored column is floored at MDL/2, so its "correlation" with another
##     mostly-censored column can be driven purely by shared non-detection. A44
##     quantifies exactly this: at PREV_THRESH = 0 more than a quarter of edges rest on
##     fewer than four joint DETECTIONS; at 0.10 not one edge does; at 0.05 only 1.7 %.
##   IS IT EXPENDABLE? Not as a censoring control; yes as a 10 % number. A44
##     (prevalence_sweep_20260904.R) recommends 0.05 -- which is exactly the Herren &
##     McMahon value: +73 compounds, z 27.9 -> 31.9, isolates 34 -> 19, 99.9 % of summed
##     toxic units, edge support essentially intact. BEFORE flipping the default here,
##     re-run scripts/unit_inventory_20260903.R at 0.05: only Basin_Elbe was tested, and
##     the small units are where StARS collapsed in A37.
## @PARAM  prevalence filter = which compounds become nodes at all ------------------------------------
PREV_THRESH <- 0.10    # minimum detection frequency (scope set by NODESET_SCOPE below)
## @PARAM  absolute-risk-driver threshold (Inostroza et al. 2024) -------------------------------------
RQ_THRESH   <- 0.02    # absolute-risk-driver threshold, after Inostroza et al. (2024)
## @PARAM  minimum sites for a basin / section unit ---------------------------------------------------
BASIN_MIN_N <- 20      # minimum sites for a river-basin-specific network


## [vRN 05.09.2026] StARS subsampling depth -----------------------------------
## Measured, not guessed: scripts/stars_depth_20260905.R, results in
## outputs_vRN/stars_depth_20260905.csv. Three basins x rep.num {20,50,100,200}
## x seeds {42,7,13}. The question is not "20 or 50" but at what depth the
## selected lambda stops depending on WHICH subsample is drawn.
##
##   basin   rep=20            rep=50 / 100 / 200
##   Elbe    0.4865, 1116 e    identical at every depth and seed
##   Danube  0.5458,  752 e    identical at every depth and seed
##   Rhine   0.6199 or 0.6712  0.6712, 617 e -- all seeds agree from 50 up
##           (789 or 617 e, depending on the seed)
##
## Two things follow.
## 1. The A23 finding that the DANUBE loses ~26 % of its edges from 20 to 50 no
##    longer reproduces. A23 (23.08.2026) predates the A36 global node set and
##    the MDL/2 floor, so it described a different node set. On the current
##    pipeline the Danube is bit-identical at every depth. Do not repeat the A23
##    caveat next to Danube numbers; cite this sweep instead.
## 2. The instability that DOES remain is the Rhine, the SMALLEST of the three
##    (63 sites), and only at rep.num = 20, where one seed in three picks a
##    different lambda and 172 edges (22 %) with it. It is gone by 50.
##
## Set to 100, not 50, for three reasons: instability clearly scales inversely
## with n and the small units (sections, reaches, n ~ 20-50 -- where StARS
## collapsed in A37) were NOT in this sweep; the cost is ~60 s per basin fit;
## and Pedro already specified StARS = 100 for the RNA-seq branch, so the two
## branches now agree. Margin is cheap here, a wrong edge set is not.
## @DECISION  StARS subsampling depth, measured not guessed (stars_depth_20260905) --------------------
STARS_REPNUM <- 100
## [vRN 03.09.2026] node-set scope ---------------------------------------------
## "global": PREV_THRESH is applied ONCE over all sites of the analysed basins
##   and the survivors become the node pool for every spatial unit. This fixes
##   the moving bar (10 % of 153 Elbe sites = 16 detections, but 10 % of a
##   20-site Rhine section = 2), which made units covering the same water carry
##   different compound sets, made node counts non-monotone in the nesting, and
##   made Q / edge count / density incomparable across units in Table 1.
## "unit": pre-03.09.2026 behaviour, filter applied inside each unit.
## @DECISION  global node set: the filter is applied ONCE, so units stay comparable -------------------
NODESET_SCOPE   <- "global"
## @PARAM  zero-variance guard inside a unit ----------------------------------------------------------
MIN_DETECT_UNIT <- 2L  # zero-variance guard: a global node needs at least this
                       # many detections in a unit, else its column is constant

## [vRN 03.09.2026] minimum sites for one along-river SECTION. 21-site Rhine thirds break
## (the middle one returns zero edges); 30-site halves do not. The number of sections per
## basin is derived from this, capped at 3 -- see make_sections() below.
## @PARAM  min sites per river section (21-site Rhine thirds returned 0 edges) ------------------------
SECTION_MIN_N <- 30

## [vRN 27.08.2026] null models -----------------------------------------------
## Supplies null_test_multi(): ER + configuration model + Chung-Lu side by side.
## See scripts/null_models_vRN.R for the full hierarchy, incl. the data-level
## (pipeline) nulls null_data() / pipeline_null().
## @NULL  THE NULL MODELS LIVE HERE -> scripts/null_models_vRN.R --------------------------------------
## [repo] the null-model hierarchy now ships with this repository
## [repo 15.09.2026] robust repo-root location. Previously `.R_DIR` was chosen by
## `if (exists("sp")) sp else "R"`, which silently does the wrong thing when an
## unrelated variable named `sp` exists in the session. Now: use `sp` only if it
## actually points at this folder; otherwise walk up from the working directory
## until run_all.R is found; last resort is the bare "R" fallback.
.locate_repo_root <- function(start = getwd()) {
  d <- tryCatch(normalizePath(start, winslash = "/", mustWork = TRUE),
                error = function(e) NULL)
  if (is.null(d)) return(NULL)
  repeat {
    if (file.exists(file.path(d, "run_all.R"))) return(d)
    parent <- dirname(d)
    if (identical(parent, d)) return(NULL)
    d <- parent
  }
}
.R_DIR <- if (exists("sp") && is.character(sp) && file.exists(file.path(sp, "null_models_vRN.R"))) {
  sp
} else {
  .root <- .locate_repo_root()
  if (!is.null(.root) && dir.exists(file.path(.root, "R"))) file.path(.root, "R") else "R"
}
source(file.path(.R_DIR, "null_models_vRN.R"))

