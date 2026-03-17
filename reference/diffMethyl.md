# Identify differentially methylated sites between conditions

The main function for differential methylation analysis in comma.
Analogous to `DESeq2::DESeq()`, `diffMethyl()` accepts a
[`commaData`](https://carl-stone.github.io/comma/reference/commaData.md)
object, fits a statistical model to each methylation site, and returns
the same object enriched with per-site test results stored as new
columns in `rowData`.

## Usage

``` r
diffMethyl(
  object,
  formula = ~condition,
  method = c("beta_binomial", "methylkit"),
  mod_type = NULL,
  min_coverage = 5L,
  p_adjust_method = "BH",
  ...
)
```

## Arguments

- object:

  A
  [`commaData`](https://carl-stone.github.io/comma/reference/commaData.md)
  object with at least two samples in distinct conditions.

- formula:

  A one-sided formula specifying the design. The RHS variable must match
  a column in `sampleInfo(object)` (e.g., `~ condition`). Default is
  `~ condition`.

- method:

  Character string selecting the statistical backend. `"beta_binomial"`
  (default) uses a quasibinomial GLM via base R. `"methylkit"` wraps
  [`methylKit::calculateDiffMeth()`](https://rdrr.io/pkg/methylKit/man/calculateDiffMeth-methods.html),
  requiring methylKit to be installed.

- mod_type:

  Character vector or `NULL`. Modification type(s) to test (e.g.,
  `"6mA"`, `c("6mA", "5mC")`). If `NULL` (default), all modification
  types present in `object` are tested.

- min_coverage:

  Integer. Minimum per-sample read depth required to include a site in
  testing. Sites where any sample has coverage below this threshold are
  treated as `NA` in that sample. Default `5L`.

- p_adjust_method:

  Character string. Multiple testing correction method, passed to
  [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html). Default `"BH"`
  (Benjamini-Hochberg). Other options: `"bonferroni"`, `"holm"`, `"BY"`,
  `"none"`.

- ...:

  Additional arguments (reserved for future use).

## Value

The input `commaData` object with additional columns in `rowData`:
`dm_pvalue`, `dm_padj`, `dm_delta_beta`, and one
`dm_mean_beta_<condition>` column per condition level. The `metadata`
slot is updated with analysis parameters and result column names.

## Details

**Statistical model (`method = "beta_binomial"`):** A per-site
quasibinomial GLM is fitted using
[`glm`](https://rdrr.io/r/stats/glm.html): \$\$
\text{Binomial}(n\_{\text{mod}},\\ n\_{\text{total}}) \sim
\text{condition} \$\$ where \\n\_{\text{mod}} = \text{round}(\beta
\times \text{coverage})\\ and \\n\_{\text{total}} = \text{coverage}\\.
The quasibinomial family accounts for overdispersion. P-values are
extracted from the Wald t-test on the contrast coefficient. This method
requires no additional packages beyond base R.

**Alternative model (`method = "methylkit"`):** Wraps
[`methylKit::calculateDiffMeth()`](https://rdrr.io/pkg/methylKit/man/calculateDiffMeth-methods.html),
which uses logistic regression. Requires methylKit to be installed
(`BiocManager::install("methylKit")`). Returns results in the same
standardised format.

**Multiple mod types:** When `mod_type = NULL` (default), all
modification types present in the object are tested independently and
results are combined. Sites of a mod type that is not being tested
receive `NA` in all `dm_*` columns.

**Small-sample note:** Differential methylation testing with very few
replicates (e.g., n = 1 per group) is mathematically possible but has
extremely low statistical power. Treat such results as exploratory only.

**Result columns added to `rowData`:**

- `dm_pvalue`:

  Raw p-value from the GLM Wald test.

- `dm_padj`:

  Adjusted p-value (Benjamini-Hochberg by default).

- `dm_delta_beta`:

  Effect size: mean methylation in the treatment group minus mean
  methylation in the reference (control) group. Positive values indicate
  hypermethylation in treatment.

- `dm_mean_beta_<condition>`:

  One column per condition level (named after the actual condition
  values), containing the per-group mean beta value for each site.

Analysis parameters and result column names are stored in
`metadata(object)$diffMethyl_params` and
`metadata(object)$diffMethyl_result_cols`.

## See also

[`results`](https://carl-stone.github.io/comma/reference/results.md) to
extract the test results as a tidy `data.frame`;
[`filterResults`](https://carl-stone.github.io/comma/reference/filterResults.md)
to filter by significance thresholds.

## Examples

``` r
data(comma_example_data)
# Test for differential 6mA methylation between conditions
dm <- diffMethyl(comma_example_data, formula = ~ condition, mod_type = "6mA")

# How many sites are significant?
rd <- as.data.frame(SummarizedExperiment::rowData(dm))
sum(rd$dm_padj < 0.05 & abs(rd$dm_delta_beta) >= 0.2, na.rm = TRUE)
#> [1] 7
```
