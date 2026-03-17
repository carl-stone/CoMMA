# Apply multiple testing correction to a vector of p-values

A thin wrapper around
[`p.adjust`](https://rdrr.io/r/stats/p.adjust.html) that operates on a
named numeric vector of p-values and returns adjusted values in the same
order. Used internally by
[`diffMethyl`](https://carl-stone.github.io/comma/reference/diffMethyl.md).

## Usage

``` r
.applyMultipleTesting(pvalues, method = "BH")
```

## Arguments

- pvalues:

  Named numeric vector of raw p-values. `NA` values are passed through
  (i.e., sites that could not be tested remain `NA` after adjustment).

- method:

  Character string. Correction method passed to
  [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html). Default `"BH"`
  (Benjamini-Hochberg FDR control). Other options include
  `"bonferroni"`, `"holm"`, `"BY"`, `"fdr"`, and `"none"`.

## Value

Named numeric vector of adjusted p-values, same length and names as
`pvalues`.
