# Metagene plot of methylation across genomic features

Computes average methylation beta values at normalized positions across
genomic features (e.g., genes from TSS to TTS) and plots the smoothed
profile. Useful for assessing whether methylation is enriched at
particular positions within a feature class.

## Usage

``` r
plot_metagene(
  object,
  feature = "gene",
  mod_type = NULL,
  motif = NULL,
  mod_context = NULL,
  n_bins = 50L
)
```

## Arguments

- object:

  A
  [`commaData`](https://carl-stone.github.io/comma/reference/commaData.md)
  object with a non-empty `annotation` slot or a user-supplied
  `features` GRanges.

- feature:

  Character string specifying the feature type to use as the reference.
  Must match a value in the `feature_type` metadata column of the
  annotation. Default `"gene"`.

- mod_type:

  Character string specifying a single modification type (e.g., `"6mA"`,
  `"5mC"`). If `NULL` (default), all modification types are used.

- motif:

  Character vector or `NULL`. If provided, only sites with matching
  sequence context motif(s) are included (e.g., `"GATC"`). If `NULL`
  (default), all motifs are included.

- n_bins:

  Positive integer. Number of equal-width bins to divide the normalized
  feature position \\\[0, 1\]\\ into. Default `50`.

## Value

A [`ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html)
object. The x-axis shows normalized position within the feature (0 =
TSS, 0.5 = midpoint, 1 = TTS); the y-axis shows mean methylation (beta).
One line is drawn per sample, colored by sample name. Dashed vertical
lines mark the TSS (0) and TTS (1).

## Details

Internally calls
[`annotateSites`](https://carl-stone.github.io/comma/reference/annotateSites.md)`(type = "metagene")`
to compute normalized positions (0 = TSS, 1 = TTS) for each methylation
site that overlaps a feature of the requested type. Sites that do not
overlap any feature are excluded from the plot. The mean beta value is
then computed within each position bin for each sample.

## See also

[`annotateSites`](https://carl-stone.github.io/comma/reference/annotateSites.md),
[`plot_genome_track`](https://carl-stone.github.io/comma/reference/plot_genome_track.md)

## Examples

``` r
data(comma_example_data)
plot_metagene(comma_example_data, feature = "gene")


# Only 6mA sites
plot_metagene(comma_example_data, feature = "gene", mod_type = "6mA")

```
