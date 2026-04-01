# Accessor for sequence context motifs present in a commaData object

Returns the sorted unique motif strings stored in
`rowData(object)$motif`. `NA` values (sites from Dorado or Megalodon
callers where motif context is unavailable) are excluded from the
result.

## Usage

``` r
motifs(object)

# S4 method for class 'commaData'
motifs(object)
```

## Arguments

- object:

  A `commaData` object.

## Value

A sorted character vector of unique non-`NA` motif strings (e.g.,
`c("CCWGG", "GATC")`). Returns `character(0)` if all motif values are
`NA` (e.g., Dorado-only data).

## Examples

``` r
data(comma_example_data)
motifs(comma_example_data)
#> [1] "CCWGG" "GATC" 
```
