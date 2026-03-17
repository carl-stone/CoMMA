# Accessor for genome size information

Returns the chromosome sizes stored in a
[`commaData`](https://carl-stone.github.io/comma/reference/commaData.md)
object.

## Usage

``` r
# S4 method for class 'commaData'
genome(x)
```

## Arguments

- x:

  A `commaData` object.

## Value

A named integer vector of chromosome sizes (chromosome name → length in
bp), or `NULL` if no genome information was provided at construction.

## Examples

``` r
data(comma_example_data)
genome(comma_example_data)
#> chr_sim 
#>  100000 
```
