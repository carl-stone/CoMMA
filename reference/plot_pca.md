# PCA of methylation profiles

Performs principal component analysis (PCA) on per-sample methylation
profiles and plots PC1 vs PC2. Useful for sample-level QC, detecting
outliers, and assessing whether biological conditions separate in
methylation space.

## Usage

``` r
plot_pca(object, mod_type = NULL, color_by = "condition", shape_by = NULL)
```

## Arguments

- object:

  A
  [`commaData`](https://carl-stone.github.io/comma/reference/commaData.md)
  object.

- mod_type:

  Character string specifying a single modification type (e.g., `"6mA"`,
  `"5mC"`). If `NULL` (default), all sites from all modification types
  are used.

- color_by:

  Character string naming a column in `sampleInfo(object)` to use for
  point color. Default `"condition"`.

- shape_by:

  Character string naming a column in `sampleInfo(object)` to use for
  point shape. If `NULL` (default), all points use the same shape.

## Value

A [`ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html)
object. PC1 and PC2 are shown on the x- and y-axes, respectively, with
percentage of variance explained shown in the axis labels. Each point
represents one sample and is labeled with its `sample_name`. Points are
colored by `color_by`.

## Details

Sites with any `NA` beta values across samples are removed before PCA to
ensure a complete data matrix. PCA is computed via
[`stats::prcomp`](https://rdrr.io/r/stats/prcomp.html) with centering
(`center = TRUE`) and without scaling (`scale. = FALSE`). A warning is
issued if fewer than three samples are present.

## See also

[`methylomeSummary`](https://carl-stone.github.io/comma/reference/methylomeSummary.md),
[`plot_methylation_distribution`](https://carl-stone.github.io/comma/reference/plot_methylation_distribution.md)

## Examples

``` r
data(comma_example_data)
plot_pca(comma_example_data)


# Color by condition, shape by replicate
plot_pca(comma_example_data, color_by = "condition")


# Only 6mA sites
plot_pca(comma_example_data, mod_type = "6mA")

```
