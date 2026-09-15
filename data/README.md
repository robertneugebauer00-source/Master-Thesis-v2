# Data

**No data files are tracked in this repository.** The `.gitignore` excludes them
deliberately -- they are third-party datasets with their own terms.

## What the pipeline expects

| File | Source |
|---|---|
| `Finckh_Carmona_2023_PANGAEA_R1.xlsx` | Finckh et al. (2024), Environ. Int. 183:108371; PANGAEA dataset [doi:10.1594/PANGAEA.960272](https://doi.pangaea.de/10.1594/PANGAEA.960272) -- see `data-raw/fetch_pangaea.R` |
| `TRIDENT_prediction_results_EC10_GRO.csv` | TRIDENT EC10 predictions |
| `MDL_CECs_summary.csv` | method detection limits, one row per compound |
| `Danube_TRIDENT_mapping.csv` | compound to SMILES map |

`data-raw/fetch_pangaea.R` downloads the PANGAEA dataset by DOI; check the sheet
layout against what `R/01_data.R` expects (`B_parameters_info`, `A_B_results`)
and point `SGH_PROJ` at the result.

## Pointing the code at your copy

`R/00_setup.R` reads the data root from an environment variable:

```r
Sys.setenv(SGH_PROJ = "C:/path/to/your/R")   # expects scripts/, raw_data/, outputs/ beneath it
```

If `SGH_PROJ` is unset it falls back to the original author path, which will not
exist on another machine -- `run_all.R` checks for the concentration file first and
tells you so before anything else runs.

Missing an optional input is not fatal: without the MDL table the code falls back to
a flat 1 ng/L non-detect floor, and without EC10 values the toxic-unit branch is
skipped with a warning.
