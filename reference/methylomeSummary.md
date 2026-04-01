# Summarize per-sample methylation and coverage distributions

Computes per-sample summary statistics for methylation beta values and
sequencing coverage in a
[`commaData`](https://carl-stone.github.io/comma/reference/commaData.md)
object. Returns a tidy `data.frame` suitable for direct use with ggplot2
or for tabular reporting.

## Usage

``` r
methylomeSummary(object, mod_type = NULL, motif = NULL)
```

## Arguments

- object:

  A
  [`commaData`](https://carl-stone.github.io/comma/reference/commaData.md)
  object.

- mod_type:

  Character string or `NULL`. If provided, only sites of the specified
  modification type (e.g., `"6mA"`) are included in the summary. If
  `NULL` (default), all modification types are summarized together.

- motif:

  Character vector or `NULL`. If provided, only sites with matching
  sequence context motif(s) are included (e.g., `"GATC"`). If `NULL`
  (default), all motifs are included.

## Value

A `data.frame` with one row per sample, containing:

- `sample_name`:

  Sample identifier.

- `condition`:

  Experimental condition, from `sampleInfo(object)$condition`.

- `mod_type`:

  The modification type summarized (`"all"` if `mod_type = NULL`).

- `n_sites`:

  Total number of sites considered.

- `n_covered`:

  Number of sites with non-`NA` methylation in this sample (i.e., sites
  above the coverage threshold).

- `mean_beta`:

  Mean beta value across covered sites.

- `median_beta`:

  Median beta value across covered sites.

- `sd_beta`:

  Standard deviation of beta values across covered sites.

- `frac_methylated`:

  Fraction of covered sites with \\\beta \> 0.5\\ (broadly methylated).

- `mean_coverage`:

  Mean sequencing depth across all sites (including sites below the
  `min_coverage` threshold, which have coverage stored as 0 or their raw
  depth).

- `median_coverage`:

  Median sequencing depth.

## See also

[`methylation`](https://carl-stone.github.io/comma/reference/methylation.md),
`coverage`,
[`sampleInfo`](https://carl-stone.github.io/comma/reference/sampleInfo.md)

## Examples

``` r
data(comma_example_data)
ms <- methylomeSummary(comma_example_data)
ms
#>   sample_name condition mod_type n_sites n_covered mean_beta median_beta
#> 1      ctrl_1   control      all     300       300 0.8678843   0.8881436
#> 2      ctrl_2   control      all     300       300 0.8728354   0.8951648
#> 3      ctrl_3   control      all     300       300 0.8781476   0.8966108
#> 4     treat_1 treatment      all     300       300 0.8135452   0.8829561
#> 5     treat_2 treatment      all     300       300 0.8136529   0.8867238
#> 6     treat_3 treatment      all     300       300 0.8004998   0.8694701
#>      sd_beta frac_methylated mean_coverage median_coverage
#> 1 0.10150073       0.9900000      76.23333            75.5
#> 2 0.09697667       0.9966667      86.88667            93.5
#> 3 0.08514535       1.0000000      81.54000            82.5
#> 4 0.19181379       0.9166667      82.13667            84.0
#> 5 0.19951744       0.9000000      78.47333            80.0
#> 6 0.21427601       0.9066667      79.93000            77.0

# Summarize only 6mA sites
ms_6mA <- methylomeSummary(comma_example_data, mod_type = "6mA")
#> Warning: missing package slot (comma) in object of class “commaData” (package info added)
ms_6mA[, c("sample_name", "condition", "mean_beta", "n_covered")]
#>   sample_name condition mean_beta n_covered
#> 1      ctrl_1   control 0.8986871       200
#> 2      ctrl_2   control 0.9002143       200
#> 3      ctrl_3   control 0.9090365       200
#> 4     treat_1 treatment 0.8189668       200
#> 5     treat_2 treatment 0.8187088       200
#> 6     treat_3 treatment 0.8001237       200
```
