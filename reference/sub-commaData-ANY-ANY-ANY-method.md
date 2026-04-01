# Subset a commaData object by sites and/or samples

Standard bracket-based subsetting. Rows correspond to methylation sites;
columns correspond to samples. The resulting object is a valid
`commaData` with all assays, rowData, colData, and custom slots updated
consistently.

## Usage

``` r
# S4 method for class 'commaData,ANY,ANY,ANY'
x[i, j, ..., drop = FALSE]
```

## Arguments

- x:

  A `commaData` object.

- i:

  Row (site) index: integer, logical, or character vector.

- j:

  Column (sample) index: integer, logical, or character vector.

- ...:

  Not used.

- drop:

  Ignored (required by generic).

## Value

A `commaData` object with the selected sites and samples.

## Examples

``` r
data(comma_example_data)
# First 50 sites, all samples
sub <- comma_example_data[1:50, ]
dim(sub)
#> [1] 50  6
```
