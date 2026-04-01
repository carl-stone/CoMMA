# Extract differential methylation results as a tidy data frame

Retrieves the per-site differential methylation statistics added by
[`diffMethyl`](https://carl-stone.github.io/comma/reference/diffMethyl.md)
and returns them as a tidy `data.frame` suitable for downstream analysis
and plotting.

## Usage

``` r
results(object, ...)

# S4 method for class 'commaData'
results(object, mod_type = NULL, motif = NULL, mod_context = NULL, ...)
```

## Arguments

- object:

  A
  [`commaData`](https://carl-stone.github.io/comma/reference/commaData.md)
  object on which
  [`diffMethyl`](https://carl-stone.github.io/comma/reference/diffMethyl.md)
  has been run.

- ...:

  Ignored (for S4 generic compatibility).

- mod_type:

  Character string or `NULL`. If provided, only sites of the specified
  modification type are returned. If `NULL` (default), results for all
  modification types are returned.

- motif:

  Character vector or `NULL`. If provided, only sites with matching
  sequence context motif(s) are returned. If `NULL` (default), all
  motifs are returned.

- mod_context:

  Character vector or `NULL`. If provided, only sites with a matching
  modification context (e.g., `"6mA_GATC"`) are returned. Use
  [`modContexts`](https://carl-stone.github.io/comma/reference/modContexts.md)`(object)`
  to see available values. Applied in addition to any `mod_type` or
  `motif` filters.

## Value

A `data.frame` with one row per methylation site, containing:

- `chrom`:

  Chromosome name.

- `position`:

  1-based genomic position.

- `strand`:

  Strand (`"+"` or `"-"`).

- `mod_type`:

  Modification type (e.g., `"6mA"`).

- `dm_pvalue`:

  Raw p-value from the statistical test.

- `dm_padj`:

  Adjusted p-value (Benjamini-Hochberg by default).

- `dm_delta_beta`:

  Effect size: mean methylation in the treatment group minus mean
  methylation in the reference group.

- `dm_mean_beta_<condition>`:

  One column per condition level with per-group mean beta values.

Any other annotation columns present in `rowData(object)` (e.g., from
[`annotateSites`](https://carl-stone.github.io/comma/reference/annotateSites.md))
are also included.

## See also

[`diffMethyl`](https://carl-stone.github.io/comma/reference/diffMethyl.md),
[`filterResults`](https://carl-stone.github.io/comma/reference/filterResults.md)

## Examples

``` r
data(comma_example_data)
dm <- diffMethyl(comma_example_data, formula = ~ condition, mod_type = "6mA")
res <- results(dm)
head(res[order(res$dm_padj), ])
#>                            chrom position strand mod_type motif mod_context
#> chr_sim:8907:-:6mA:GATC  chr_sim     8907      -      6mA  GATC    6mA_GATC
#> chr_sim:52014:+:6mA:GATC chr_sim    52014      +      6mA  GATC    6mA_GATC
#> chr_sim:69527:+:6mA:GATC chr_sim    69527      +      6mA  GATC    6mA_GATC
#> chr_sim:72824:-:6mA:GATC chr_sim    72824      -      6mA  GATC    6mA_GATC
#> chr_sim:62293:-:6mA:GATC chr_sim    62293      -      6mA  GATC    6mA_GATC
#> chr_sim:9028:-:6mA:GATC  chr_sim     9028      -      6mA  GATC    6mA_GATC
#>                          is_diff    dm_pvalue     dm_padj dm_delta_beta
#> chr_sim:8907:-:6mA:GATC     TRUE 1.209684e-04 0.009737468    -0.6097199
#> chr_sim:52014:+:6mA:GATC    TRUE 6.241567e-05 0.009737468    -0.7591100
#> chr_sim:69527:+:6mA:GATC    TRUE 1.640873e-04 0.009737468    -0.6702704
#> chr_sim:72824:-:6mA:GATC   FALSE 1.947494e-04 0.009737468    -0.1353157
#> chr_sim:62293:-:6mA:GATC    TRUE 3.217300e-04 0.012869199    -0.6522237
#> chr_sim:9028:-:6mA:GATC     TRUE 5.508523e-04 0.018361742    -0.6850134
#>                          dm_mean_beta_control dm_mean_beta_treatment
#> chr_sim:8907:-:6mA:GATC             0.9362931              0.3265732
#> chr_sim:52014:+:6mA:GATC            0.9284467              0.1693367
#> chr_sim:69527:+:6mA:GATC            0.9604202              0.2901498
#> chr_sim:72824:-:6mA:GATC            0.9809048              0.8455891
#> chr_sim:62293:-:6mA:GATC            0.9627737              0.3105501
#> chr_sim:9028:-:6mA:GATC             0.9345386              0.2495252
```
