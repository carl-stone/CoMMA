# Extract differential methylation results as a tidy data frame

Retrieves the per-site differential methylation statistics added by
[`diffMethyl`](https://carl-stone.github.io/comma/reference/diffMethyl.md)
and returns them as a tidy `data.frame` suitable for downstream analysis
and plotting.

## Usage

``` r
results(object, ...)

# S4 method for class 'commaData'
results(object, mod_type = NULL, ...)
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
#>                       chrom position strand mod_type is_diff    dm_pvalue
#> chr_sim:8907:-:6mA  chr_sim     8907      -      6mA    TRUE 1.209684e-04
#> chr_sim:52014:+:6mA chr_sim    52014      +      6mA    TRUE 6.241567e-05
#> chr_sim:69527:+:6mA chr_sim    69527      +      6mA    TRUE 1.640873e-04
#> chr_sim:72824:-:6mA chr_sim    72824      -      6mA   FALSE 1.947494e-04
#> chr_sim:62293:-:6mA chr_sim    62293      -      6mA    TRUE 3.217300e-04
#> chr_sim:9028:-:6mA  chr_sim     9028      -      6mA    TRUE 5.508523e-04
#>                         dm_padj dm_delta_beta dm_mean_beta_control
#> chr_sim:8907:-:6mA  0.009737468    -0.6097199            0.9362931
#> chr_sim:52014:+:6mA 0.009737468    -0.7591100            0.9284467
#> chr_sim:69527:+:6mA 0.009737468    -0.6702704            0.9604202
#> chr_sim:72824:-:6mA 0.009737468    -0.1353157            0.9809048
#> chr_sim:62293:-:6mA 0.012869199    -0.6522237            0.9627737
#> chr_sim:9028:-:6mA  0.018361742    -0.6850134            0.9345386
#>                     dm_mean_beta_treatment
#> chr_sim:8907:-:6mA               0.3265732
#> chr_sim:52014:+:6mA              0.1693367
#> chr_sim:69527:+:6mA              0.2901498
#> chr_sim:72824:-:6mA              0.8455891
#> chr_sim:62293:-:6mA              0.3105501
#> chr_sim:9028:-:6mA               0.2495252
```
