# Subset a commaData object by condition, modification type, or chromosome

A convenience function for filtering a
[`commaData`](https://carl-stone.github.io/comma/reference/commaData.md)
object by common criteria. For arbitrary index-based subsetting, use
`[`.

## Usage

``` r
subset(x, ...)

# S4 method for class 'commaData'
subset(x, mod_type = NULL, condition = NULL, chrom = NULL, motif = NULL, ...)
```

## Arguments

- x:

  A `commaData` object.

- ...:

  Ignored.

- mod_type:

  Character vector or `NULL`. If provided, only sites with a matching
  modification type are kept (e.g., `"6mA"`).

- condition:

  Character vector or `NULL`. If provided, only samples matching the
  specified condition(s) are kept.

- chrom:

  Character vector or `NULL`. If provided, only sites on the specified
  chromosome(s) are kept.

- motif:

  Character vector or `NULL`. If provided, only sites with a matching
  sequence context motif are kept (e.g., `"GATC"`). Sites with `NA`
  motif values are excluded when this filter is active. Use
  [`motifs`](https://carl-stone.github.io/comma/reference/motifs.md) to
  see which motifs are present.

## Value

A `commaData` object containing only the selected sites and samples.

## Examples

``` r
data(comma_example_data)
# Only 6mA sites
six_ma <- subset(comma_example_data, mod_type = "6mA")
#> Warning: missing package slot (comma) in object of class “commaData” (package info added)
modTypes(six_ma)
#> [1] "6mA"

# Only GATC-context sites
gatc <- subset(comma_example_data, motif = "GATC")
#> Warning: missing package slot (comma) in object of class “commaData” (package info added)
nrow(gatc)
#> [1] 200
```
