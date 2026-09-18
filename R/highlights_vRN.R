## =============================================================================
## highlights_vRN.R -- a navigable index of the important parts of the workflow
## Robert Neugebauer. Created 14.09.2026.
##
## The tags are plain R comments, so nothing here changes any result. The point
## is that RStudio can jump straight to them, and the index can be printed or
## exported for a meeting.
##
## TAG VOCABULARY -- put one of these anywhere in a comment line:
##   #@KEY       the load-bearing step of a section
##   #@NULL      anything to do with the null hypothesis / H0
##   #@PARAM     a tunable parameter, and the justification for its value
##   #@DECISION  a choice that was made, and why (provenance, references)
##   #@CAVEAT    a limitation that must be stated when the result is reported
##   #@TODO      open item
##
## USAGE
##   source("scripts/highlights_vRN.R")
##   hl()                  # every tag, grouped
##   hl("NULL")            # only the null-hypothesis tags
##   hl_jump("NULL", 2)    # open the Rmd at the 2nd @NULL tag
##   hl_sections()         # headings + chunk labels, i.e. the map of the file
##   hl_md()               # export the index as markdown
## =============================================================================

## [repo 15.09.2026] resolve against the repo ROOT (found by walking up to
## run_all.R), not against whatever the working directory happens to be --
## the old relative paths silently broke outside the repo root.
.hl_root <- function(start = getwd()) {
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
.HL_ROOT <- .hl_root()
HL_FILE <- Sys.getenv("SGH_HL_FILE",
                      unset = if (!is.null(.HL_ROOT)) file.path(.HL_ROOT, "analysis", "03_mainline.Rmd") else
                              file.path("analysis", "03_mainline.Rmd"))
HL_TAGS <- c("KEY", "NULL", "PARAM", "DECISION", "CAVEAT", "TODO")
HL_OUT  <- if (!is.null(.HL_ROOT)) file.path(.HL_ROOT, "outputs", "script_highlights.md") else
           file.path("outputs", "script_highlights.md")

## ---------------------------------------------------------------------------
## scan
## ---------------------------------------------------------------------------
hl_scan <- function(file = HL_FILE, tags = HL_TAGS) {
  x <- readLines(file, warn = FALSE)
  hits <- grep("#+ *@[A-Z]", x)                  # matches  #@KEY  and  ## @KEY
  if (!length(hits)) return(data.frame(line = integer(0), tag = character(0),
                                       text = character(0), stringsAsFactors = FALSE))
  res <- do.call(rbind, lapply(hits, function(i) {
    s    <- x[i]
    p    <- regexpr("@[A-Z]", s)
    rest <- substring(s, p + 1)
    tag  <- toupper(sub(" .*", "", rest))
    txt  <- trimws(sub("[-=]{3,}$", "", trimws(substring(rest, nchar(tag) + 1))))
    data.frame(line = i, tag = tag, text = txt, stringsAsFactors = FALSE)
  }))
  res <- res[res$tag %in% tags, , drop = FALSE]
  res[order(match(res$tag, tags), res$line), , drop = FALSE]
}

## ---------------------------------------------------------------------------
## print
## ---------------------------------------------------------------------------
hl <- function(tag = NULL, file = HL_FILE) {
  d <- hl_scan(file)
  if (!is.null(tag)) d <- d[d$tag == toupper(tag), , drop = FALSE]
  if (!nrow(d)) { message("no tags found -- add #@KEY / #@NULL / ... to the script"); return(invisible(d)) }
  for (tg in unique(d$tag)) {
    dd <- d[d$tag == tg, , drop = FALSE]
    cat("\n== @", tg, "  (", nrow(dd), ")\n", sep = "")
    for (k in seq_len(nrow(dd)))
      cat(sprintf("   L%-5d %s\n", dd$line[k], substr(dd$text[k], 1, 96)))
  }
  cat("\n")
  invisible(d)
}

## ---------------------------------------------------------------------------
## jump (RStudio only)
## ---------------------------------------------------------------------------
hl_jump <- function(tag, n = 1, file = HL_FILE) {
  d <- hl_scan(file)
  d <- d[d$tag == toupper(tag), , drop = FALSE]
  if (!nrow(d))    stop("no @", toupper(tag), " tags in ", basename(file))
  if (n > nrow(d)) stop("only ", nrow(d), " @", toupper(tag), " tags")
  ## [repo 15.09.2026] guard: degrade gracefully outside RStudio
  if (!requireNamespace("rstudioapi", quietly = TRUE)) {
    message("rstudioapi not available -- the tag is at ", file, ":", d$line[n])
    return(invisible(d[n, ]))
  }
  rstudioapi::navigateToFile(file, line = d$line[n])
  invisible(d[n, ])
}

## ---------------------------------------------------------------------------
## the map of the file: markdown headings + code chunk labels
## ---------------------------------------------------------------------------
hl_sections <- function(file = HL_FILE, print = TRUE) {
  x  <- readLines(file, warn = FALSE)
  tx <- trimws(x)
  fence    <- startsWith(tx, "```")
  in_chunk <- (cumsum(fence) - fence) %% 2 == 1   # TRUE for lines inside a code chunk
  is_head  <- grepl("^#{1,3} ", tx) & !in_chunk & !fence
  is_chunk <- startsWith(tx, "```{r")
  idx <- sort(c(which(is_head), which(is_chunk)))
  d <- data.frame(line = idx,
                  kind = ifelse(is_chunk[idx], "chunk", "heading"),
                  text = tx[idx], stringsAsFactors = FALSE)
  if (print) for (k in seq_len(nrow(d)))
    cat(sprintf("L%-5d %-8s %s\n", d$line[k], d$kind[k], substr(d$text[k], 1, 92)))
  invisible(d)
}

## ---------------------------------------------------------------------------
## export the index (markdown, opens anywhere)
## ---------------------------------------------------------------------------
hl_md <- function(file = HL_FILE, out = HL_OUT) {
  d <- hl_scan(file)
  dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
  ln <- c(paste("# Script highlights --", basename(file)),
          paste("Generated", format(Sys.time(), "%d.%m.%Y %H:%M")), "")
  for (tg in unique(d$tag)) {
    dd <- d[d$tag == tg, , drop = FALSE]
    ln <- c(ln, paste0("## @", tg), "")
    ln <- c(ln, paste0("- **L", dd$line, "** ", dd$text), "")
  }
  writeLines(ln, out)
  message("written: ", out)
  invisible(out)
}

message("highlights_vRN.R loaded -- try hl(), hl_sections(), hl_jump(\"NULL\")")
