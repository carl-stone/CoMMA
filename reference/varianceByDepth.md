# Methylation variance as a function of sequencing depth

For each coverage level (or coverage bin), computes the variance of
methylation beta values across sites at that depth. This is useful for
diagnosing whether low-coverage sites have inflated methylation variance
and for setting appropriate coverage thresholds.

## Usage

``` r
varianceByDepth(object, coverage_bins = NULL, mod_type = NULL)
```

## Arguments

- object:

  A
  [`commaData`](https://carl-stone.github.io/comma/reference/commaData.md)
  object.

- coverage_bins:

  Integer vector specifying the coverage levels to include. If `NULL`
  (default), all unique coverage levels observed across all samples are
  used. Useful to pass `5:30` to focus on a specific depth range.

- mod_type:

  Character string or `NULL`. If provided, only sites of the specified
  modification type are included. Default: `NULL` (all types).

## Value

A `data.frame` with one row per (coverage level, sample), containing:

- `coverage`:

  Sequencing depth (integer).

- `sample_name`:

  Sample identifier.

- `variance`:

  Variance of beta values at sites with exactly this coverage level.
  `NA` if fewer than 2 sites are at this level.

- `n_sites`:

  Number of sites at this coverage level.

## See also

[`coverageDepth`](https://carl-stone.github.io/comma/reference/coverageDepth.md),
[`methylomeSummary`](https://carl-stone.github.io/comma/reference/methylomeSummary.md)

## Examples

``` r
data(comma_example_data)
vd <- varianceByDepth(comma_example_data, coverage_bins = 5:30)
head(vd)
#>   coverage sample_name    variance n_sites
#> 1        5      ctrl_1          NA       0
#> 2        6      ctrl_1          NA       0
#> 3        7      ctrl_1          NA       0
#> 4        8      ctrl_1          NA       0
#> 5        9      ctrl_1          NA       0
#> 6       10      ctrl_1 0.002581064       2
```
