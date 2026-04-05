# Windowed sequencing depth across the genome

Bins the genome into non-overlapping windows and computes the average
(or median) sequencing depth in each window for each sample. Returns a
tidy `data.frame` suitable for plotting or downstream QC analysis.

## Usage

``` r
coverageDepth(
  object,
  window,
  method = c("mean", "median"),
  log2_transform = FALSE
)
```

## Arguments

- object:

  A
  [`commaData`](https://carl-stone.github.io/comma/reference/commaData.md)
  object.

- window:

  Positive integer. Window size in base pairs.

- method:

  Character string. Aggregation method within each window. One of
  `"mean"` (default) or `"median"`.

- log2_transform:

  Logical. If `TRUE`, the depth values are log2-transformed (using
  \\log2(depth + 1)\\ to handle zeros). Default: `FALSE`.

## Value

A `data.frame` with one row per (chromosome window, sample), containing:

- `chrom`:

  Chromosome name.

- `window_start`:

  First base of the window (1-based).

- `window_end`:

  Last base of the window (1-based).

- `sample_name`:

  Sample identifier.

- `depth`:

  Mean or median sequencing depth in the window.

- `log2_depth`:

  Log2-transformed depth (only present if `log2_transform = TRUE`).

## Details

Depth is computed only at positions with observed methylation sites.
Windows with no sites have `depth = NA`.

If genome size information is stored in `genome(object)`, windows are
sized to fit the chromosomes exactly (the last window may be smaller
than `window`). If genome information is absent, only the range spanned
by observed sites is covered.

## See also

[`varianceByDepth`](https://carl-stone.github.io/comma/reference/varianceByDepth.md),
[`methylomeSummary`](https://carl-stone.github.io/comma/reference/methylomeSummary.md)

## Examples

``` r
data(comma_example_data)
cd <- coverageDepth(comma_example_data, window = 10000L)
head(cd)
#>     chrom window_start window_end sample_name    depth
#> 1 chr_sim            1      10000      ctrl_1 78.98276
#> 2 chr_sim        10001      20000      ctrl_1 88.43396
#> 3 chr_sim        20001      30000      ctrl_1 82.53333
#> 4 chr_sim        30001      40000      ctrl_1 70.52632
#> 5 chr_sim        40001      50000      ctrl_1 71.08929
#> 6 chr_sim        50001      60000      ctrl_1 72.08696
```
