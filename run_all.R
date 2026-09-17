## =============================================================================
## run_all.R -- run the whole pipeline, in order, from one place.
## Robert Neugebauer, MSc thesis. Created 14.09.2026.
##
## Works in BOTH layouts without editing:
##   repo layout   R/00_setup.R ... analysis/03_mainline.Rmd
##   flat layout   00_setup.R ... 03_mainline.Rmd all in one folder
##
## Run it with the working directory set to the repo root (or the split folder):
##   setwd("C:/Users/rober/git/sgh-chemical-networks"); source("run_all.R")
##
## Nothing here changes any analysis. It only decides WHAT runs and in WHICH order,
## times each step, and writes a log so a failed run says where it stopped.
## =============================================================================

## ---- 1. what to run --------------------------------------------------------
## Set a step to FALSE to skip it. The first three are the shared foundation and
## every analysis step needs them, so leave them TRUE unless they are already
## loaded in the session. `lioness` (05_lioness.R) is also side-effect-free
## function definitions -- only S6 calls them, but sourcing costs nothing.
STEPS <- c(
  setup     = TRUE,   # 00_setup.R      palette, packages, paths, ALL parameters
  data      = TRUE,   # 01_data.R       PANGAEA, MDL/2 floor, toxic units, units
  functions = TRUE,   # 02_functions.R  prep / fit / metrics / modules / risk / plot
  lioness   = TRUE,   # 05_lioness.R    LIONESS/BONOBO single-sample network functions
  mainline  = TRUE,   # 03_mainline.Rmd THE PIPELINE -- Tables 1-4, Figures 1-4
  S1        = FALSE,  # method justification   (slow: prevalence sweep, netCompare)
  S2        = FALSE,  # SBM / DC-SBM           (slow: 326 lines of model fitting)
  S3        = FALSE,  # cross-basin rarefaction (slow: R refits per basin)
  S4        = FALSE,  # null-model hierarchy   (slow: N draws x 10 null families)
  S5        = FALSE,  # the SGH test itself    (stress index, Q ~ stress)
  S6        = FALSE,  # single-sample networks (LIONESS/BONOBO per-site nets)
  S7        = FALSE,  # site-level tests (pooled inference, UDF proxy, nulls; needs S6 output)
  S8        = FALSE,  # backbone sensitivity (needs S6 run with LIONESS_GLASSO <- TRUE)
  S9        = FALSE,  # paired chemical+microbial simulation, known ground truth (slow: power grid)
  utils     = FALSE   # 99_utils.R      figure refresh + HTML report builder
)

## "render" knits each .Rmd to HTML (the report you would show someone).
## "purl"   extracts the code and sources it -- same results, no HTML, much faster.
MODE          <- "render"
STOP_ON_ERROR <- TRUE      # FALSE = carry on and report which steps failed

## ---- 2. locate everything --------------------------------------------------
R_DIR  <- if (dir.exists("R"))        "R"        else "."
AN_DIR <- if (dir.exists("analysis")) "analysis" else "."
OUT    <- if (dir.exists("outputs"))  "outputs"  else "."
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
LOG <- file.path(OUT, paste0("run_all_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))

FILES <- list(
  setup     = file.path(R_DIR,  "00_setup.R"),
  data      = file.path(R_DIR,  "01_data.R"),
  functions = file.path(R_DIR,  "02_functions.R"),
  lioness   = file.path(R_DIR,  "05_lioness.R"),
  mainline  = file.path(AN_DIR, "03_mainline.Rmd"),
  S1        = file.path(AN_DIR, "S1_method_justification.Rmd"),
  S2        = file.path(AN_DIR, "S2_alternative_communities.Rmd"),
  S3        = file.path(AN_DIR, "S3_crossbasin.Rmd"),
  S4        = file.path(AN_DIR, "S4_null_hierarchy.Rmd"),
  S5        = file.path(AN_DIR, "S5_stress_gradient.Rmd"),
  S6        = file.path(AN_DIR, "S6_single_networks.Rmd"),
  S7        = file.path(AN_DIR, "S7_site_level_tests.Rmd"),
  S8        = file.path(AN_DIR, "S8_backbone_sensitivity.Rmd"),
  S9        = file.path(AN_DIR, "S9_simulation.Rmd"),
  utils     = file.path(R_DIR,  "99_utils.R")
)

## ---- 3. helpers ------------------------------------------------------------
say <- function(...) { msg <- paste0(...); cat(msg, "\n", sep = ""); cat(msg, "\n", sep = "", file = LOG, append = TRUE) }

hms <- function(sec) sprintf("%02d:%02d", as.integer(sec) %/% 60, as.integer(sec) %% 60)

run_one <- function(key) {
  path <- FILES[[key]]
  if (!file.exists(path)) return(list(status = "MISSING", secs = 0, msg = path))
  t0 <- Sys.time()
  res <- tryCatch({
    if (grepl("[.]Rmd$", path)) {
      if (MODE == "render") {
        if (!requireNamespace("rmarkdown", quietly = TRUE)) stop("rmarkdown not installed; set MODE to purl")
        rmarkdown::render(path, output_dir = OUT, envir = globalenv(), quiet = TRUE)
      } else {
        tmp <- tempfile(fileext = ".R")
        knitr::purl(path, output = tmp, quiet = TRUE)
        source(tmp, echo = FALSE)
      }
    } else {
      source(path, echo = FALSE)
      ## [repo 15.09.2026] 99_utils.R is side-effect-free on source (defines
      ## functions only); drive its entry point explicitly.
      if (key == "utils" && exists("main", envir = globalenv(), inherits = FALSE) &&
          is.function(get("main", envir = globalenv())))
        get("main", envir = globalenv())()
    }
    "OK"
  }, error = function(e) paste("ERROR:", conditionMessage(e)))
  list(status = res, secs = as.numeric(difftime(Sys.time(), t0, units = "secs")), msg = path)
}

## ---- 4. pre-flight ---------------------------------------------------------
say("=== run_all.R  ", format(Sys.time(), "%d.%m.%Y %H:%M:%S"), " ===")
say("working directory : ", getwd())
say("layout            : ", if (R_DIR == "R") "repo (R/ + analysis/)" else "flat")
say("mode              : ", MODE)
say("steps requested   : ", paste(names(STEPS)[STEPS], collapse = ", "))

miss <- names(STEPS)[STEPS & !file.exists(unlist(FILES[names(STEPS)]))]
if (length(miss)) say("!! FILES NOT FOUND for: ", paste(miss, collapse = ", "))

## the data step is the one that fails for a reason worth naming
if (isTRUE(STEPS[["data"]])) {
  .root <- Sys.getenv("SGH_PROJ", unset = "C:/Users/rober/OneDrive/Masterkram/R")
  .xl   <- file.path(.root, "scripts", "Finckh_Carmona_2023_PANGAEA_R1.xlsx")
  if (!file.exists(.xl))
    say("!! concentration file not found: ", .xl,
        "   -- point the SGH_PROJ environment variable at your data tree")
}

## ---- 5. run ----------------------------------------------------------------
results <- list()
for (key in names(STEPS)) {
  if (!isTRUE(STEPS[[key]])) { say("-- skip  ", key); next }
  say("-> ", key, "  (", basename(FILES[[key]]), ")")
  r <- run_one(key)
  results[[key]] <- r
  say("   ", r$status, "   ", hms(r$secs))
  if (STOP_ON_ERROR && grepl("^ERROR", r$status)) { say("!! stopping (STOP_ON_ERROR = TRUE)"); break }
}

## ---- 6. summary ------------------------------------------------------------
if (length(results)) {
  summ <- data.frame(step   = names(results),
                     file   = basename(unlist(lapply(results, `[[`, "msg"))),
                     status = unlist(lapply(results, `[[`, "status")),
                     mmss   = hms(unlist(lapply(results, `[[`, "secs"))),
                     stringsAsFactors = FALSE, row.names = NULL)
  say("")
  say("=== summary ===")
  for (i in seq_len(nrow(summ))) say(sprintf("  %-10s %-32s %-8s %s", summ$step[i], summ$file[i], substr(summ$status[i], 1, 8), summ$mmss[i]))
  say("total: ", hms(sum(unlist(lapply(results, `[[`, "secs")))))
  say("log:   ", normalizePath(LOG, winslash = "/", mustWork = FALSE))
  write.csv(summ, file.path(OUT, "run_all_summary.csv"), row.names = FALSE)
  print(summ, row.names = FALSE)
} else say("nothing ran -- every step was set to FALSE")

invisible(results)
