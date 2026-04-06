# Plot coverage depth distribution

Produces a histogram of sequencing depth (coverage) across sites for
each sample in a
[`commaData`](https://carl-stone.github.io/comma/reference/commaData.md)
object. Useful for QC to assess whether coverage is sufficient and
consistent across samples.

## Usage

``` r
plot_coverage(
  object,
  mod_type = NULL,
  motif = NULL,
  mod_context = NULL,
  per_sample = TRUE
)
```

## Arguments

- object:

  A
  [`commaData`](https://carl-stone.github.io/comma/reference/commaData.md)
  object.

- mod_type:

  Character string specifying a single modification type (e.g., `"6mA"`,
  `"5mC"`). If `NULL` (default), all sites from all modification types
  are included.

- motif:

  Character vector or `NULL`. If provided, only sites with matching
  sequence context motif(s) are included (e.g., `"GATC"`). If `NULL`
  (default), all motifs are included.

- mod_context:

  Character vector or `NULL`. If provided, only sites whose
  `mod_context` rowData column matches one of the supplied values are
  included. If `NULL` (default), all modification contexts are included.

- per_sample:

  Logical. If `TRUE` (default), the plot is faceted by sample, producing
  one histogram panel per sample. If `FALSE`, all samples are overlaid
  on a single plot with per-sample colors.

## Value

A [`ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html)
object. The x-axis shows coverage depth on a log10 scale; the y-axis
shows the number of sites at each depth. A vertical dashed line marks
the median coverage per sample (when `per_sample = TRUE`) or across all
samples (when `per_sample = FALSE`). Sites with `NA` coverage are
silently excluded.

## See also

[`coverageDepth`](https://carl-stone.github.io/comma/reference/coverageDepth.md),
[`varianceByDepth`](https://carl-stone.github.io/comma/reference/varianceByDepth.md),
[`plot_methylation_distribution`](https://carl-stone.github.io/comma/reference/plot_methylation_distribution.md)

## Examples

``` r
data(comma_example_data)
plot_coverage(comma_example_data)


# Overlay all samples on one plot
plot_coverage(comma_example_data, per_sample = FALSE)


# One modification type only
plot_coverage(comma_example_data, mod_type = "6mA")

```
