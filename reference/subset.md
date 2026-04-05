# Subset a commaData object by condition, modification type, or chromosome

A convenience function for filtering a
[`commaData`](https://carl-stone.github.io/comma/reference/commaData.md)
object by common criteria. For arbitrary index-based subsetting, use
`[`.

## Usage

``` r
subset(x, ...)

# S4 method for class 'commaData'
subset(
  x,
  mod_type = NULL,
  condition = NULL,
  chrom = NULL,
  motif = NULL,
  mod_context = NULL,
  ...
)
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

- mod_context:

  Character vector or `NULL`. If provided, only sites with a matching
  modification context are kept (e.g., `"6mA_GATC"`, `"5mC_CCWGG"`). A
  `mod_context` value is `paste(mod_type, motif, sep = "_")` when motif
  is available, or just `mod_type` for Dorado/Megalodon data. Use
  [`modContexts`](https://carl-stone.github.io/comma/reference/modContexts.md)
  to see which contexts are present. When provided, this filter is
  applied in addition to (ANDed with) any `mod_type` or `motif` filters.

## Value

A `commaData` object containing only the selected sites and samples.

## Examples

``` r
data(comma_example_data)
# Only 6mA sites
six_ma <- subset(comma_example_data, mod_type = "6mA")
modTypes(six_ma)
#> [1] "6mA"

# Only GATC-context sites
gatc <- subset(comma_example_data, motif = "GATC")
nrow(gatc)
#> [1] 393

# Filter by mod_context (equivalent to the above for modkit data)
gatc2 <- subset(comma_example_data, mod_context = "6mA_GATC")
nrow(gatc2)
#> [1] 393
```
