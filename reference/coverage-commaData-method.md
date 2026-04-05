# Accessor for the sequencing coverage (read depth) matrix

Retrieves the sites × samples matrix of read depth from a
[`commaData`](https://carl-stone.github.io/comma/reference/commaData.md)
object.

## Usage

``` r
# S4 method for class 'commaData'
coverage(x, shift = 0L, width = NULL, weight = 1L, ...)
```

## Arguments

- x:

  A `commaData` object.

- shift:

  Not used; inherited from the
  [`IRanges::coverage`](https://rdrr.io/pkg/IRanges/man/coverage-methods.html)
  generic.

- width:

  Not used; inherited from the
  [`IRanges::coverage`](https://rdrr.io/pkg/IRanges/man/coverage-methods.html)
  generic.

- weight:

  Not used; inherited from the
  [`IRanges::coverage`](https://rdrr.io/pkg/IRanges/man/coverage-methods.html)
  generic.

- ...:

  Not used.

## Value

An integer matrix with rows corresponding to methylation sites and
columns corresponding to samples.

## See also

[`methylation`](https://carl-stone.github.io/comma/reference/methylation.md),
[`siteInfo`](https://carl-stone.github.io/comma/reference/siteInfo.md)

## Examples

``` r
data(comma_example_data)
cov <- coverage(comma_example_data)
summary(as.vector(cov))
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>   10.00   44.00   79.00   79.31  114.00  150.00 
```
