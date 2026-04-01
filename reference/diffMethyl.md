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
  method = c("beta_binomial", "methylkit", "limma", "quasi_f"),
  mod_type = NULL,
  motif = NULL,
  min_coverage = 5L,
  alpha = 0.5,
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
  (default) uses a quasibinomial GLM via base R. `"quasi_f"` applies
  empirical Bayes shrinkage of quasibinomial dispersions via
  [`squeezeVar`](https://rdrr.io/pkg/limma/man/squeezeVar.html)
  (quasi-likelihood F-test; count-data EB, recommended for small n).
  Requires limma. `"limma"` applies empirical Bayes variance shrinkage
  via [`eBayes`](https://rdrr.io/pkg/limma/man/ebayes.html) on
  M-value-transformed data; recommended when replicates are few (n \< 3
  per group). Requires limma. `"methylkit"` wraps
  [`methylKit::calculateDiffMeth()`](https://rdrr.io/pkg/methylKit/man/calculateDiffMeth-methods.html),
  requiring methylKit to be installed.

- mod_type:

  Character vector or `NULL`. Modification type(s) to test (e.g.,
  `"6mA"`, `c("6mA", "5mC")`). If `NULL` (default), all modification
  types present in `object` are tested.

- motif:

  Character vector or `NULL`. If provided, only sites with matching
  sequence context motif(s) are tested (e.g., `"GATC"`). Uses
  [`motifs`](https://carl-stone.github.io/comma/reference/motifs.md) to
  validate the requested values. If `NULL` (default), all motifs
  (including `NA`) are included.

- min_coverage:

  Integer. Minimum per-sample read depth required to include a site in
  testing. Sites where any sample has coverage below this threshold are
  treated as `NA` in that sample. Default `5L`.

- alpha:

  Positive numeric pseudocount used to compute M-values when
  `method = "limma"`: \\M = \log_2((n\_{\mathrm{mod}} + \alpha) /
  (n\_{\mathrm{unmod}} + \alpha))\\. Default `0.5` (a Beta(0.5, 0.5)
  prior). Ignored for other methods.

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

**Alternative model (`method = "quasi_f"`):** A two-pass extension of
`"beta_binomial"` that adds empirical Bayes shrinkage of the per-site
quasibinomial dispersion estimates, analogous to the quasi-likelihood
F-test of edgeR. Pass 1 fits the same quasibinomial GLM per site and
collects the per-site dispersion \\\hat\phi_j\\ and residual df
\\df_j\\. Pass 2 calls
[`squeezeVar`](https://rdrr.io/pkg/limma/man/squeezeVar.html) to
estimate a log-normal prior on \\\\\hat\phi_j\\\\ and compute posterior
estimates \\\\\tilde\phi_j\\\\. Pass 3 recomputes the t-statistic using
\\\tilde\phi_j\\ and evaluates it against a t-distribution with \\d_0 +
df_j\\ degrees of freedom (where \\d_0\\ is the estimated prior df).
Requires limma.

**Alternative model (`method = "limma"`):** Beta values are transformed
to M-values via \\M = \log_2((n\_{\mathrm{mod}} + \alpha) /
(n\_{\mathrm{unmod}} + \alpha))\\, then
[`lmFit`](https://rdrr.io/pkg/limma/man/lmFit.html) fits an OLS model
per site and [`eBayes`](https://rdrr.io/pkg/limma/man/ebayes.html)
applies empirical Bayes variance shrinkage — borrowing information
across all sites to stabilize the per-site variance estimate. This gives
substantially more power than `"beta_binomial"` when replicates are few
(n \< 3 per group). Requires limma (`BiocManager::install("limma")`).
Effect sizes are reported on the original beta scale.

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
