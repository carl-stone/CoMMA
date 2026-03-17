# Per-site quasibinomial GLM for differential methylation

Fits a per-site quasibinomial generalized linear model to test for
differential methylation between conditions. This function is an
internal engine called by
[`diffMethyl`](https://carl-stone.github.io/comma/reference/diffMethyl.md)
when `method = "beta_binomial"`.

## Usage

``` r
.betaBinomialTest(methyl_mat, coverage_mat, coldata, formula)
```

## Arguments

- methyl_mat:

  Numeric matrix (sites × samples) of beta values. `NA` indicates
  below-coverage sites.

- coverage_mat:

  Integer matrix (sites × samples) of read depths.

- coldata:

  `data.frame` with at least one column matching the RHS variable in
  `formula` (typically `condition`).

- formula:

  One-sided formula specifying the design (e.g., `~ condition`).

## Value

A `data.frame` with one row per site (same row order as `methyl_mat`),
containing:

- `pvalue`:

  Raw p-value from the Wald test on the condition coefficient. `NA` for
  untestable sites.

- `delta_beta`:

  Estimated effect size: mean beta in the treatment group minus mean
  beta in the reference (control) group. Uses the model matrix reference
  level as the baseline. `NA` for untestable sites.

- `mean_beta_<level>`:

  One column per unique condition level, named `mean_beta_<level>`
  (e.g., `mean_beta_control`, `mean_beta_treatment`). Contains the
  per-group observed mean beta value computed directly from non-NA data.

## Details

The model treats methylation counts as overdispersed binomial data: for
each site, the number of modified reads is modelled as
\\\text{Binomial}(n\_{\text{coverage}}, p)\\ with a quasibinomial family
to account for overdispersion. The formula is fitted site-by-site using
[`glm`](https://rdrr.io/r/stats/glm.html).

**Small-sample caveat:** With few replicates (e.g., 2 control + 1
treatment), residual degrees of freedom will be very small (df = 1 in
this case). The quasibinomial dispersion estimate and the resulting
p-values remain mathematically valid but will have very low power.
Interpret results from n \< 3 per group with caution.

Sites are skipped (result is `NA`) if:

- All samples have `NA` beta values (low coverage).

- Fewer than 2 non-NA samples are present.

- The GLM fails to converge or the contrast coefficient cannot be
  extracted (e.g., singular fit, only one condition level present after
  NA removal).
