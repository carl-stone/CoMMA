# CoMMA development news

## Unreleased

- Refactored `writeBED()` to a file-path driven API with explicit input arguments, input validation, and deterministic BED9 output (removes hard-coded local paths and global `WT_average` reliance).
- Improved roxygen documentation for exported annotation and rolling-summary functions to document explicit input/return schemas.
- Replaced placeholder test coverage with unit tests for exported annotation/rolling/BED APIs and key internal differential-methylation helper behavior.
- Declared missing runtime dependencies used by package code (`tibble`, `zoo`) in `DESCRIPTION`.
- Extended canonical site-table validation to support base-agnostic methylation workflows: `mod_base` is now backfilled with a configurable default (`"6mA"`) for legacy inputs and `motif` is supported as optional metadata with default `"GATC"`.
- Updated `run_differential_methylation()` to add explicit `mod_base` and `motif` filters (defaults preserve historical 6mA/GATC behavior) while keeping the output result schema stable.
- Updated `run_differential_methylation()` with a new `context_label` parameter (default `"CpG"`) that is passed through to methylKit `methylRaw(..., context = ...)` metadata for explicit bacterial modification/motif context labeling.
- Added exported `run_differential_methylation()` high-level wrapper for differential methylation using methylKit with manuscript-default coverage filtering, median normalization, SLIM correction, and default significance thresholds.
