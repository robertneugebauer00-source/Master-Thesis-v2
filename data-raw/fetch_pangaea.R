## =============================================================================
## data-raw/fetch_pangaea.R -- fetch the concentration dataset from PANGAEA
## [repo 15.09.2026]
##
## Downloads the published dataset behind
## Finckh_Carmona_2023_PANGAEA_R1.xlsx, the file the pipeline expects
## (R/00_setup.R: file.path(SGH_PROJ, "scripts", "Finckh_Carmona_2023_PANGAEA_R1.xlsx"),
## sheets "B_parameters_info" and "A_B_results").
##
## CITATION (always quote the citation above the data when using them):
##   Finckh, Sven; Carmona, Eric; et al. (2024): Concentrations of 504 organic
##   micropollutants in 445 European streams by LC-HRMS target screening
##   (2016-2019). PANGAEA, https://doi.org/10.1594/PANGAEA.960272
##   Published as: Finckh, S. et al. (2024): Mapping chemical footprints of
##   organic micropollutants in European streams. Environment International,
##   183:108371, https://doi.org/10.1016/j.envint.2023.108371
##   -- PANGAEA etiquette: "Always quote citation above when using data."
##
## RUN:  Rscript data-raw/fetch_pangaea.R      (from the repo root)
## =============================================================================

## ---- locate the repo root (walk up from this script until run_all.R) ---------
.args     <- commandArgs(trailingOnly = FALSE)
.file_arg <- sub("^--file=", "", grep("^--file=", .args, value = TRUE))
.start    <- if (length(.file_arg)) dirname(normalizePath(.file_arg[1], winslash = "/")) else
             getwd()
ROOT <- .start
repeat {
  if (file.exists(file.path(ROOT, "run_all.R"))) break
  .parent <- dirname(ROOT)
  if (identical(.parent, ROOT)) stop("fetch_pangaea.R: could not locate repo root above ", .start)
  ROOT <- .parent
}

raw_dir <- file.path(ROOT, "data-raw")
out_dir <- file.path(raw_dir, "pangaea_960272")
zipfile <- file.path(raw_dir, "pangaea_960272.zip")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## ---- (a) download the dataset bundle -----------------------------------------
url <- "https://doi.pangaea.de/10.1594/PANGAEA.960272?format=zip"
cat("downloading", url, "\n ->", zipfile, "\n")
download.file(url, zipfile, mode = "wb")   # mode "wb" is required for binary/zip on Windows

## ---- (b) unzip ----------------------------------------------------------------
cat("unzipping ->", out_dir, "\n")
unzip(zipfile, exdir = out_dir)

## ---- (c) show what arrived ------------------------------------------------------
arrived <- list.files(out_dir, recursive = TRUE)
cat("archive contents:\n")
cat(paste0("  - ", arrived), sep = "\n")

## ---- (d) what to do next ---------------------------------------------------------
cat("
================================================================================
NEXT STEPS -- the download is NOT yet plug-compatible with the pipeline.

1. The pipeline expects ONE xlsx workbook,
     <SGH_PROJ>/scripts/Finckh_Carmona_2023_PANGAEA_R1.xlsx
   with exactly these two sheets (see R/00_setup.R and R/01_data.R):
     - 'B_parameters_info'  (site metadata; header row contains 'Sample code',
                             'Latitude', 'Longitude', 'Country', 'River basin')
     - 'A_B_results'        (compounds in rows, sites in columns; a header row
                             whose first cell is 'Name' and whose site columns
                             are named EUS_*)
   PANGAEA ships tab-delimited text, one table per subset. VERIFY the sheet
   names and layouts match the above; if they do not, export/rename the sheets
   accordingly (e.g. read the tabs, write a two-sheet xlsx with openxlsx).

2. Point the pipeline at the result by setting SGH_PROJ to the data root that
   holds scripts/Finckh_Carmona_2023_PANGAEA_R1.xlsx, e.g.
     Sys.setenv(SGH_PROJ = \"<path-to-your-data-root>\")   # see data/README.md

3. The remaining inputs are NOT on PANGAEA and must be obtained separately
   (see data/README.md):
     - MDL_CECs_summary.csv                 (method detection limits)
     - TRIDENT_prediction_results_EC10_GRO.csv  and  Danube_TRIDENT_mapping.csv
       (EC10 toxicity predictions + compound-to-SMILES map, for the TU branch)
   Without them the pipeline still runs: it falls back to a flat 1 ng/L
   non-detect floor and skips the toxic-unit branch with a warning.

Citation reminder: always quote the PANGAEA citation (header of this script)
when using these data.
================================================================================
")

invisible(arrived)
