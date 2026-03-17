# Run differential methylation via methylKit

An internal wrapper that uses methylKit to test for differential
methylation. Called by
[`diffMethyl`](https://carl-stone.github.io/comma/reference/diffMethyl.md)
when `method = "methylkit"`.

## Usage

``` r
.runMethylKit(methyl_mat, coverage_mat, coldata, formula)
```

## Arguments

- methyl_mat:

  Numeric matrix (sites × samples) of beta values.

- coverage_mat:

  Integer matrix (sites × samples) of read depths.

- coldata:

  `data.frame` with at least one column matching the RHS variable in
  `formula`.

- formula:

  One-sided formula (e.g., `~ condition`).

## Value

A `data.frame` with the same columns as
[`.betaBinomialTest()`](https://carl-stone.github.io/comma/reference/dot-betaBinomialTest.md):
`pvalue`, `delta_beta`, and one `mean_beta_<level>` column per condition
level. Row names are site keys.

## Details

methylKit must be installed (it is listed in `Suggests`). If it is not
available, this function stops with an informative message.

The wrapper converts the methylation and coverage matrices from a
[`commaData`](https://carl-stone.github.io/comma/reference/commaData.md)
object into the format expected by
[`methylKit::methylRawList`](https://rdrr.io/pkg/methylKit/man/methylRawList-class.html),
runs
[`methylKit::unite()`](https://rdrr.io/pkg/methylKit/man/unite-methods.html)
and
[`methylKit::calculateDiffMeth()`](https://rdrr.io/pkg/methylKit/man/calculateDiffMeth-methods.html),
and returns results in the same standardised format as
[`.betaBinomialTest()`](https://carl-stone.github.io/comma/reference/dot-betaBinomialTest.md).

Only the first RHS variable of `formula` is used as the grouping
variable. Complex formulas with interactions or batch terms are not
currently supported by this wrapper.
