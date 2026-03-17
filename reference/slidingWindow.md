# Sliding window methylation summary along the genome

Computes a per-position sliding window statistic (median or mean) of
methylation beta values for each sample in a
[`commaData`](https://carl-stone.github.io/comma/reference/commaData.md)
object. The genome size for each chromosome is read from
`genome(object)`, so no organism-specific values are ever hardcoded.
Circular genome wrap-around is supported via
[`.circularIndex()`](https://carl-stone.github.io/comma/reference/dot-circularIndex.md).

## Usage

``` r
slidingWindow(
  object,
  window,
  stat = c("median", "mean"),
  mod_type = NULL,
  circular = TRUE
)
```

## Arguments

- object:

  A
  [`commaData`](https://carl-stone.github.io/comma/reference/commaData.md)
  object. Must have genome size information in `genome(object)` (i.e.,
  it must have been constructed with a `genome` argument).

- window:

  Positive integer. Window size in base pairs. The smoothed value at
  position \\p\\ is computed from positions \\\[p - \lfloor w/2
  \rfloor,\\ p + \lfloor w/2 \rfloor\]\\.

- stat:

  Character string. Summary statistic to apply within each window. One
  of `"median"` (default) or `"mean"`.

- mod_type:

  Character string or `NULL`. If provided, only sites of the specified
  modification type (e.g., `"6mA"`) are included in the smoothing. If
  `NULL` (default), all sites are used.

- circular:

  Logical. If `TRUE` (default), positions at the ends of each chromosome
  are wrapped around so that the window at position 1 can draw from
  positions near the chromosome end, and vice versa. Appropriate for
  circular bacterial chromosomes.

## Value

A `data.frame` with one row per (chromosome, position, sample)
combination, containing:

- `chrom`:

  Chromosome name (character).

- `position`:

  Genomic position, 1-based (integer).

- `sample_name`:

  Sample identifier (character).

- `window_stat`:

  Smoothed beta value (numeric). The column is named `window_median` or
  `window_mean` depending on `stat`.

## Details

Because most positions in a genome have no methylation site, the beta
value vector for each chromosome-sample pair is sparse (mostly `NA`).
[`zoo::rollapply`](https://rdrr.io/pkg/zoo/man/rollapply.html) is called
with `na.rm = TRUE` so that windows spanning regions with no data still
produce a value where at least one site is present.

Positions where every site in the window is `NA` (i.e., no coverage)
remain `NA` in the output.

## See also

[`methylation`](https://carl-stone.github.io/comma/reference/methylation.md),
`genome`,
[`methylomeSummary`](https://carl-stone.github.io/comma/reference/methylomeSummary.md)

## Examples

``` r
data(comma_example_data)
# \donttest{
sw <- slidingWindow(comma_example_data, window = 5000L)
head(sw)
#>     chrom position sample_name window_median
#> 1 chr_sim        1      ctrl_1     0.8581977
#> 2 chr_sim        2      ctrl_1     0.8581977
#> 3 chr_sim        3      ctrl_1     0.8581977
#> 4 chr_sim        4      ctrl_1     0.8581977
#> 5 chr_sim        5      ctrl_1     0.8581977
#> 6 chr_sim        6      ctrl_1     0.8581977
# Filter to one sample
sw_ctrl1 <- sw[sw$sample_name == "ctrl_1", ]
# }
```
