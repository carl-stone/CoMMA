# Plot methylation beta value distributions

Produces a density plot of methylation beta values (0–1) for each sample
in a
[`commaData`](https://carl-stone.github.io/comma/reference/commaData.md)
object. Useful for QC and for comparing methylation level distributions
across samples and modification types.

## Usage

``` r
plot_methylation_distribution(
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

  Character string specifying a single modification type to plot (e.g.,
  `"6mA"`, `"5mC"`). If `NULL` (default), all modification types are
  included and the plot is faceted by `mod_type`.

- motif:

  Character vector or `NULL`. If provided, only sites with matching
  sequence context motif(s) are included (e.g., `"GATC"`). If `NULL`
  (default), all motifs are included.

- per_sample:

  Logical. If `TRUE` (default), a separate density curve is drawn for
  each sample. If `FALSE`, a single aggregate density curve is drawn per
  modification type.

## Value

A [`ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html)
object. The x-axis shows beta values (0 = unmethylated, 1 = fully
methylated); the y-axis shows kernel density. When `per_sample = TRUE`,
curves are colored by `sample_name`. When multiple modification types
are present (and `mod_type = NULL`), the plot is faceted by `mod_type`.
Sites with `NA` beta values (below coverage threshold) are silently
excluded.

## See also

[`methylomeSummary`](https://carl-stone.github.io/comma/reference/methylomeSummary.md),
[`plot_coverage`](https://carl-stone.github.io/comma/reference/plot_coverage.md)

## Examples

``` r
data(comma_example_data)
plot_methylation_distribution(comma_example_data)


# One modification type only
plot_methylation_distribution(comma_example_data, mod_type = "6mA")


# Aggregate across samples
plot_methylation_distribution(comma_example_data, per_sample = FALSE)

```
