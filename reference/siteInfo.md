# Accessor for per-site metadata

Returns the per-site metadata table from a
[`commaData`](https://carl-stone.github.io/comma/reference/commaData.md)
object. Equivalent to `rowData(object)` but returns a plain
`data.frame`.

## Usage

``` r
siteInfo(object)

# S4 method for class 'commaData'
siteInfo(object)
```

## Arguments

- object:

  A `commaData` object.

## Value

A `data.frame` with one row per methylation site. Always contains
columns `chrom`, `position`, `strand`, and `mod_type`. May contain
additional annotation columns added by
[`annotateSites()`](https://carl-stone.github.io/comma/reference/annotateSites.md).

## See also

[`methylation`](https://carl-stone.github.io/comma/reference/methylation.md),
[`modTypes`](https://carl-stone.github.io/comma/reference/modTypes.md)

## Examples

``` r
data(comma_example_data)
head(siteInfo(comma_example_data))
#> DataFrame with 6 rows and 5 columns
#>                          chrom  position      strand    mod_type   is_diff
#>                    <character> <integer> <character> <character> <logical>
#> chr_sim:444:+:6mA      chr_sim       444           +         6mA      TRUE
#> chr_sim:1072:+:6mA     chr_sim      1072           +         6mA     FALSE
#> chr_sim:1600:-:6mA     chr_sim      1600           -         6mA     FALSE
#> chr_sim:2176:+:6mA     chr_sim      2176           +         6mA      TRUE
#> chr_sim:3565:+:6mA     chr_sim      3565           +         6mA     FALSE
#> chr_sim:4767:-:6mA     chr_sim      4767           -         6mA     FALSE
```
