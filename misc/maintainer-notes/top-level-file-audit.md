# Top-level file audit (2026-02-23)

This note records the disposition of files that previously lived at repository root.

| Original root file | Classification | Action |
|---|---|---|
| `.DS_Store` | legacy analysis artifact | Removed (macOS metadata file; not source). |
| `WT_6mA_Mg.txt` | reproducible example fixture | Moved to `inst/extdata/WT_6mA_Mg.txt`. |
| `WT_6mA_all_callers.txt` | reproducible example fixture | Moved to `inst/extdata/WT_6mA_all_callers.txt`. |
| `all_site_annotations.txt` | reproducible example fixture | Moved to `inst/extdata/all_site_annotations.txt`. |
| `all_site_annotations_60p.txt` | reproducible example fixture | Moved to `inst/extdata/all_site_annotations_60p.txt`. |
| `functions.R` | legacy analysis artifact | Moved to `misc/legacy/functions.R`. |
| `methylKitGATC.R` | legacy analysis artifact | Moved to `misc/legacy/methylKitGATC.R`. |
| `testscript.R` | legacy analysis artifact | Moved to `misc/legacy/testscript.R`. |

## Package-source locations (policy reminder)

Package source code and shipping assets should be kept under canonical package paths:

- `R/`, `tests/`, `vignettes/`, `inst/` for active package code/data
- `inst/extdata/` for reproducible input fixtures
- `tests/testthat/testdata/` for test-only fixtures
- `misc/legacy/` for historical analysis scripts that are not package API
