# AGENTS.md — CoMMA (Comparison of Microbial Methylated Adenines)
**Read this file first.** It defines “done,” constraints, and how to work safely in parallel.

## Mission
CoMMA is an R package for **comparative bacterial methylomics**: ingest per-site base-modification calls (Nanopore/PacBio-derived callers), annotate sites with genomic context, and perform robust differential methylation analysis between groups/conditions/clones.

**Near-term scope:** keep the package name as-is (CoMMA), but **design features to be base-agnostic** (support 6mA / 5mC / 4mC in bacteria) without breaking existing 6mA-at-GATC workflows.

## Non-negotiables (Definition of Done)
A PR is “mergeable” only if:
1. `R CMD check --no-manual --as-cran` passes (no ERROR/WARNING; NOTES must be explained/justified).
2. `devtools::test()` passes (add/extend tests for any bugfix or new feature).
3. `devtools::document()` is run and committed (roxygen + NAMESPACE stay consistent).
4. Public API changes are **intentional**:
   - If you change an exported function signature or output schema: update docs + add tests + note it in `NEWS.md`.
5. Do not introduce new heavy dependencies casually. Prefer:
   - core deps in `Imports`
   - optional features in `Suggests` with graceful fallback.

## Core pipeline contract (do not silently change)
These behaviors are part of the intended method:
- Differential methylation uses **methylKit** as the underlying engine (logistic regression; SLIM correction).
- Sites with insufficient coverage are excluded (**<10 reads** cutoff).
- Coverage is normalized by **median coverage**.
- Differential calls use **q-value < 0.05** and **|percent methylation difference| > 10%** thresholds.
- If a motif site is mutated in any sample/clone, that site is removed from differential analysis so only shared sites are compared.

If you propose changing any of the above, do it as an explicit, reviewable design decision:
- add a parameter with a default that preserves current behavior
- document it
- add tests demonstrating the old and new behavior

## Repo hygiene rules for parallel “agent swarm”
- One PR per focused change set.
- Keep PRs small enough to review.
- Never mix refactors with new features in the same PR unless unavoidable.
- Prefer “surgical” diffs: minimal surface area, strong tests.

Recommended branch pattern:
- `rehab/mainline` is the integration branch
- feature branches: `rehab/<topic>-<shortslug>`

## What the package should DO (product goals)
### A) Ingest methylation calls (base-agnostic)
Support a standardized internal representation of per-site modification calls:
- Required columns (internal canonical schema):
  - `seqname`, `pos`, `strand` (or equivalent)
  - `mod_base` (e.g., "6mA", "5mC", "4mC")
  - `motif` (optional; e.g., "GATC")
  - `n_mod`, `n_total` (counts) OR `beta` + coverage
  - sample identifiers: `sample_id`, `group` (condition)
- Provide import helpers for common caller outputs (Megalodon/Tombo/others) as best-effort converters.
- Always preserve raw counts when available; compute `beta = n_mod / n_total` consistently.

### B) Annotate sites with genomic features
- Accept feature annotations via BED/GFF/GenBank-derived tables.
- Provide helpers to attach annotations (genic/intergenic, promoter elements, TF binding sites, repeats, prophage, IS elements, etc.).
- Annotation must not require E. coli-specific databases; those can be optional adapters.

### C) Differential methylation
- Provide a high-level function that:
  - filters low coverage
  - normalizes coverage
  - runs methylKit differential methylation
  - returns a tidy result table + a stable set of plotting-ready objects
- Allow:
  - pairwise comparisons (A vs B)
  - multi-group designs (if supported by underlying methods)
  - replicate-aware analysis

### D) Exploratory analyses (optional modules)
- Ordination / distance comparisons on normalized beta-values.
- GO enrichment as a *downstream* analysis (should be optional via `Suggests`).

## Commands agents must run
From repo root:

### Fast iteration
- `R -q -e 'devtools::document()'`
- `R -q -e 'devtools::test()'`

### Pre-PR gate
- `R -q -e 'devtools::check(args=c(\"--no-manual\",\"--as-cran\"))'`

### Optional quality gates (only if already set up in repo)
- `R -q -e 'lintr::lint_package()'`   (if lintr configured)
- `R -q -e 'pkgdown::build_site()'`   (if pkgdown configured)

## Testing expectations (be strict)
- Every exported function should have:
  - at least one unit test for normal behavior
  - tests for input validation (bad types, missing columns)
  - tests for stable column names / output schema
- Add small synthetic datasets in `tests/testthat/_snaps` or as fixtures under `tests/testthat/testdata/`.
- Avoid giant binary fixtures.

## Documentation expectations
- Roxygen docs must include:
  - expected input schema
  - return value schema
  - examples that run fast
- Add/maintain:
  - `README` quickstart
  - at least one “end-to-end” vignette demonstrating:
    ingest → annotate → differential → summarize/plot

## Work allocation (“swarm” task menu)
Pick ONE task per PR. Good default ordering:

1) **Stabilize package infrastructure**
   - Make checks reliable
   - Fix NAMESPACE/roxygen drift
   - Add minimal tests scaffold

2) **Codify the internal data schema**
   - Create a canonical “site table” class or validation helpers
   - Add converters for current formats used in the manuscript

3) **Harden differential methylation wrapper**
   - Ensure the pipeline contract is implemented exactly
   - Add tests that confirm the filtering/normalization/cutoffs

4) **Generalize beyond 6mA/GATC**
   - Extend schema to `mod_base`
   - Ensure default behavior remains 6mA/GATC-compatible
   - Add one example/test for 5mC or 4mC with a toy motif

5) **Annotation adapters**
   - Generic BED/GFF annotation support
   - Optional E. coli RegulonDB-style adapters (behind `Suggests`)

6) **Docs + vignette polish**
   - Turn the manuscript workflow into a clean vignette
   - Keep runtime short; use toy data

## Things NOT to do without explicit instruction
- Do not rename the package now.
- Do not rewrite the entire API.
- Do not introduce heavy dependencies (Bioconductor, huge genomics stacks) unless gated behind `Suggests` and justified.
- Do not change statistical defaults/cutoffs silently.

## If you are unsure
Open a PR that only adds:
- tests documenting current behavior
- small refactors that reduce complexity without changing outputs

Prefer “make it correct and testable” over “make it fancy.”
