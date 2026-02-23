# Contributing to CoMMA

## Directory conventions

Please keep repository contents in canonical package locations:

- `R/`: exported and internal package functions.
- `tests/testthat/`: unit tests.
- `tests/testthat/testdata/`: test-only fixtures.
- `vignettes/`: end-to-end package workflows.
- `inst/extdata/`: reproducible example input files shipped with the package.
- `misc/legacy/`: historical analysis scripts and non-package artifacts kept for reference.

Avoid adding ad hoc analysis files to repository root. If a file is needed for package behavior, tests, or user examples, place it in one of the directories above.
