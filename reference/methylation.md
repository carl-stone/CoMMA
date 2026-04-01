# Accessor for the methylation (beta value) matrix

Retrieves the sites × samples matrix of methylation beta values from a
[`commaData`](https://carl-stone.github.io/comma/reference/commaData.md)
object. Values are in the range 0–1 (proportion of reads called
methylated). Sites below the `min_coverage` threshold set at object
creation are `NA`.

## Usage

``` r
methylation(object)

# S4 method for class 'commaData'
methylation(object)
```

## Arguments

- object:

  A `commaData` object.

## Value

A numeric matrix with rows corresponding to methylation sites and
columns corresponding to samples. Rownames are site keys
(`"chrom:position:strand:mod_type"`); column names are sample names.

## See also

[`coverage`](https://carl-stone.github.io/comma/reference/coverage-commaData-method.md),
[`siteInfo`](https://carl-stone.github.io/comma/reference/siteInfo.md),
[`sampleInfo`](https://carl-stone.github.io/comma/reference/sampleInfo.md)

## Examples

``` r
data(comma_example_data)
m <- methylation(comma_example_data)
dim(m)
#> [1] 300   6
head(m)
#>                            ctrl_1    ctrl_2    ctrl_3    treat_1   treat_2
#> chr_sim:444:+:6mA:GATC  0.9048367 0.8237897 0.8517928 0.30677235 0.4768195
#> chr_sim:1072:+:6mA:GATC 0.8128778 0.8501235 0.9732947 0.74361241 0.8354918
#> chr_sim:1600:-:6mA:GATC 0.8581977 0.8124665 0.9528697 0.94747731 0.9132681
#> chr_sim:2176:+:6mA:GATC 0.9033107 0.8708096 0.9598866 0.02648254 0.5052399
#> chr_sim:3565:+:6mA:GATC 0.8558892 0.9512447 0.9277963 0.88531486 0.9560292
#> chr_sim:4767:-:6mA:GATC 0.7848818 0.8599388 0.9546899 0.91202669 0.9682563
#>                           treat_3
#> chr_sim:444:+:6mA:GATC  0.4124745
#> chr_sim:1072:+:6mA:GATC 0.9364048
#> chr_sim:1600:-:6mA:GATC 0.9212486
#> chr_sim:2176:+:6mA:GATC 0.2795380
#> chr_sim:3565:+:6mA:GATC 0.9166728
#> chr_sim:4767:-:6mA:GATC 0.9142556
```
