# Accessor for per-sample metadata

Returns the per-sample metadata table from a
[`commaData`](https://carl-stone.github.io/comma/reference/commaData.md)
object. Equivalent to `colData(object)` but returns a plain `data.frame`
for ease of use.

## Usage

``` r
sampleInfo(object)

# S4 method for class 'commaData'
sampleInfo(object)
```

## Arguments

- object:

  A `commaData` object.

## Value

A `data.frame` with one row per sample. Always contains columns
`sample_name`, `condition`, and `replicate`. May contain additional
columns such as `caller` and `file_path`.

## See also

[`siteInfo`](https://carl-stone.github.io/comma/reference/siteInfo.md),
[`modTypes`](https://carl-stone.github.io/comma/reference/modTypes.md)

## Examples

``` r
data(comma_example_data)
sampleInfo(comma_example_data)
#>         sample_name condition replicate caller
#> ctrl_1       ctrl_1   control         1 modkit
#> ctrl_2       ctrl_2   control         2 modkit
#> ctrl_3       ctrl_3   control         3 modkit
#> treat_1     treat_1 treatment         1 modkit
#> treat_2     treat_2 treatment         2 modkit
#> treat_3     treat_3 treatment         3 modkit
```
