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
columns `chrom`, `position`, `strand`, `mod_type`, `motif` (the sequence
context; `NA` for Dorado/Megalodon callers), and `mod_context` (the
composite modification context, e.g., `"6mA_GATC"`). May contain
additional annotation columns added by
[`annotateSites()`](https://carl-stone.github.io/comma/reference/annotateSites.md).

## See also

[`methylation`](https://carl-stone.github.io/comma/reference/methylation.md),
[`modTypes`](https://carl-stone.github.io/comma/reference/modTypes.md)

## Examples

``` r
data(comma_example_data)
head(siteInfo(comma_example_data))
#> DataFrame with 6 rows and 7 columns
#>                               chrom  position      strand    mod_type
#>                         <character> <integer> <character> <character>
#> chr_sim:443:+:6mA:GATC      chr_sim       443           +         6mA
#> chr_sim:512:+:6mA:GATC      chr_sim       512           +         6mA
#> chr_sim:1024:+:6mA:GATC     chr_sim      1024           +         6mA
#> chr_sim:1073:+:6mA:GATC     chr_sim      1073           +         6mA
#> chr_sim:1536:-:6mA:GATC     chr_sim      1536           -         6mA
#> chr_sim:1602:+:6mA:GATC     chr_sim      1602           +         6mA
#>                               motif mod_context   is_diff
#>                         <character> <character> <logical>
#> chr_sim:443:+:6mA:GATC         GATC    6mA_GATC     FALSE
#> chr_sim:512:+:6mA:GATC         GATC    6mA_GATC     FALSE
#> chr_sim:1024:+:6mA:GATC        GATC    6mA_GATC     FALSE
#> chr_sim:1073:+:6mA:GATC        GATC    6mA_GATC     FALSE
#> chr_sim:1536:-:6mA:GATC        GATC    6mA_GATC     FALSE
#> chr_sim:1602:+:6mA:GATC        GATC    6mA_GATC      TRUE
```
