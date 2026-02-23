# CoMMA development news

## Unreleased

- Extended canonical site-table validation to support base-agnostic methylation workflows: `mod_base` is now backfilled with a configurable default (`"6mA"`) for legacy inputs and `motif` is supported as optional metadata with default `"GATC"`.
- Updated `run_differential_methylation()` to add explicit `mod_base` and `motif` filters (defaults preserve historical 6mA/GATC behavior) while keeping the output result schema stable.
- Added exported `run_differential_methylation()` high-level wrapper for differential methylation using methylKit with manuscript-default coverage filtering, median normalization, SLIM correction, and default significance thresholds.
