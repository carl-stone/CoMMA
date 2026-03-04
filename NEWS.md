# CoMMA development news

## Unreleased

- **Breaking:** `writeBED()` `chrom` parameter no longer defaults to `"U00096.3"` (E. coli K-12). Callers must now supply `chrom` explicitly.
- **Breaking:** `methylRollingMedian()` `genome_size` parameter no longer defaults to `4641652` (E. coli). Callers must now supply `genome_size` explicitly.
- **Breaking:** `methylRollingMedian(..., method = "fast")` output column renamed from `mean_methyl` to `med_methyl` (values have always been medians; the name was wrong).
- `run_differential_methylation()` default `context_label` changed from `"CpG"` to `"6mA_GATC"`. This label is only metadata embedded in methylKit objects and does not affect statistical results. Pass `context_label = "CpG"` explicitly to restore the old label.
- `run_differential_methylation()` gains a `reference_group` parameter. When set, the named group is assigned `treatment = 0` (reference/denominator) in methylKit, making the sign of `percent_diff` explicit and user-controlled. Default `NULL` preserves the previous alphabetical assignment.

- Added a legacy compatibility roadmap: `annotateMethylSites()` now wraps canonical annotation helpers internally (`validate_site_table()` + `normalize_annotation_table()` + `annotate_sites_with_features()`) while preserving prior output markers/columns.
- Soft-deprecated legacy helpers `annotateTSS()`, `methylRollingMean()`, `calculateMethylSiteDepth()`, and `varByCoverage()` with once-per-session lifecycle warnings and explicit migration targets.
- Added README migration guidance mapping legacy calls to canonical workflows.
- Updated `validate_site_table()` to support two explicit measurement modes: count mode (`n_mod` + `n_total`) or fraction mode (`beta` + `coverage`), with strict mode validation and canonical output columns including backfilled integer counts in fraction mode.
- Refactored `writeBED()` to a file-path driven API with explicit input arguments, input validation, and deterministic BED9 output (removes hard-coded local paths and global `WT_average` reliance).
- Improved roxygen documentation for exported annotation and rolling-summary functions to document explicit input/return schemas.
- Replaced placeholder test coverage with unit tests for exported annotation/rolling/BED APIs and key internal differential-methylation helper behavior.
- Declared missing runtime dependencies used by package code (`tibble`, `zoo`) in `DESCRIPTION`.
- Extended canonical site-table validation to support base-agnostic methylation workflows: `mod_base` is now backfilled with a configurable default (`"6mA"`) for legacy inputs and `motif` is supported as optional metadata with default `"GATC"`.
- Updated `run_differential_methylation()` to add explicit `mod_base` and `motif` filters (defaults preserve historical 6mA/GATC behavior) while keeping the output result schema stable.
- Updated `run_differential_methylation()` with a new `context_label` parameter (default `"CpG"`) that is passed through to methylKit `methylRaw(..., context = ...)` metadata for explicit bacterial modification/motif context labeling.
- Added exported `run_differential_methylation()` high-level wrapper for differential methylation using methylKit with manuscript-default coverage filtering, median normalization, SLIM correction, and default significance thresholds.
