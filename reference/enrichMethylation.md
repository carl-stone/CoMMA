# Gene set enrichment analysis of differential methylation results

Maps per-site differential methylation statistics up to the gene level
and runs GO and/or KEGG enrichment analysis using clusterProfiler.
Supports both over-representation analysis (ORA) and gene set enrichment
analysis (GSEA). Works with model organisms via standard Bioconductor
databases, or with any organism via a custom `TERM2GENE` mapping.

## Usage

``` r
enrichMethylation(
  object,
  method = "ora",
  OrgDb = NULL,
  keyType = "SYMBOL",
  ont = "BP",
  organism = NULL,
  TERM2GENE = NULL,
  TERM2NAME = NULL,
  gene_col = "feature_names",
  feature_type = "gene",
  padj_threshold = 0.05,
  delta_beta_threshold = 0.1,
  score_metric = "combined",
  gene_score_agg = "max",
  mod_type = NULL,
  mod_context = NULL,
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  minGSSize = 10L,
  maxGSSize = 500L
)
```

## Arguments

- object:

  A
  [`commaData`](https://carl-stone.github.io/comma/reference/commaData.md)
  object on which
  [`diffMethyl`](https://carl-stone.github.io/comma/reference/diffMethyl.md)
  **and**
  [`annotateSites`](https://carl-stone.github.io/comma/reference/annotateSites.md)
  have been run.

- method:

  Character vector; one or both of `"ora"` (over- representation
  analysis) and `"gsea"` (gene set enrichment analysis). Default
  `"ora"`.

- OrgDb:

  A Bioconductor `OrgDb` annotation object for GO enrichment (e.g.,
  `org.EcK12.eg.db`). If `NULL` and no `TERM2GENE` is provided, GO
  analysis is skipped.

- keyType:

  Character string specifying the key type used for gene IDs in `OrgDb`
  (e.g., `"SYMBOL"`, `"ENTREZID"`). Must match the identifiers stored in
  `gene_col`. Default `"SYMBOL"`.

- ont:

  Character string specifying the GO ontology to use when `OrgDb` is
  provided. One of `"BP"`, `"MF"`, `"CC"`, or `"ALL"`. Default `"BP"`.

- organism:

  Character string; KEGG organism code (e.g., `"eco"`). If `NULL`, KEGG
  analysis is skipped. Requires internet access.

- TERM2GENE:

  A two-column `data.frame` with columns `term` and `gene` mapping
  pathway/term identifiers to gene identifiers. When provided, this is
  used for GO analysis via
  [`enricher`](https://rdrr.io/pkg/clusterProfiler/man/enricher.html)
  (ORA) or [`GSEA`](https://rdrr.io/pkg/clusterProfiler/man/GSEA.html)
  and takes precedence over `OrgDb`. Ignored for KEGG (use `organism`
  for KEGG).

- TERM2NAME:

  An optional two-column `data.frame` with columns `term` and `name`
  providing human-readable descriptions of term IDs in `TERM2GENE`.

- gene_col:

  Character string; the `rowData` column containing gene identifiers per
  site (a `CharacterList` or `list` column added by
  [`annotateSites`](https://carl-stone.github.io/comma/reference/annotateSites.md)).
  Default `"feature_names"`.

- feature_type:

  Character vector or `NULL`. When non-`NULL`, only sites whose
  `feature_types` annotation (added by
  [`annotateSites`](https://carl-stone.github.io/comma/reference/annotateSites.md))
  contains at least one entry matching `feature_type` are used, and
  `gene_col` is subset to the matching names only. Default `"gene"` —
  restricts enrichment to gene-body sites, which is required by most
  pathway databases. Set to `NULL` to include all annotated sites
  regardless of feature type.

- padj_threshold:

  Numeric; adjusted p-value threshold for classifying a site as
  differentially methylated in ORA. Default `0.05`.

- delta_beta_threshold:

  Numeric; minimum absolute effect size (\\\|\Delta\beta\|\\) for
  classifying a site as differentially methylated in ORA. Default `0.1`.

- score_metric:

  Character string controlling the per-site scoring function used for
  GSEA ranking. One of `"combined"` (default), `"padj"`, or
  `"delta_beta"`. See Details.

- gene_score_agg:

  Character string; how to aggregate per-site scores to a per-gene score
  when multiple sites map to the same gene. Either `"max"` (default;
  largest absolute score, sign preserved) or `"mean"`.

- mod_type:

  Character string or `NULL`; passed to
  [`results`](https://carl-stone.github.io/comma/reference/results.md)
  for modification-type filtering before enrichment.

- mod_context:

  Character string or `NULL`; passed to
  [`results`](https://carl-stone.github.io/comma/reference/results.md)
  for modification-context filtering.

- pvalueCutoff:

  Numeric; p-value cutoff passed to clusterProfiler. Default `0.05`.

- qvalueCutoff:

  Numeric; q-value cutoff passed to clusterProfiler (ORA only). Default
  `0.2`.

- minGSSize:

  Integer; minimum gene set size passed to clusterProfiler. Default
  `10`.

- maxGSSize:

  Integer; maximum gene set size passed to clusterProfiler. Default
  `500`.

## Value

A named list with elements `$go` and `$kegg` (either or both may be
`NULL` if the corresponding analysis was not requested or not possible).
When a single `method` is requested, each element is a clusterProfiler
`enrichResult` (ORA) or `gseaResult` (GSEA), or `NULL`. When both
`"ora"` and `"gsea"` are requested, each element is itself a named list
with `$ora` and `$gsea` slots.

## Prerequisites

Before calling `enrichMethylation()`, you must:

1.  Run
    [`diffMethyl`](https://carl-stone.github.io/comma/reference/diffMethyl.md)
    to compute per-site `dm_padj` and `dm_delta_beta` values.

2.  Run
    [`annotateSites`](https://carl-stone.github.io/comma/reference/annotateSites.md)
    (with `type = "overlap"`) to assign gene identifiers to each site.
    The resulting `feature_names` column (a `CharacterList`) is used by
    default.

## Gene-to-pathway mapping

Supply at least one of the following:

- `TERM2GENE`:

  A two-column `data.frame` with columns `term` (pathway/GO term ID) and
  `gene` (gene identifier matching values in `gene_col`). Used for the
  GO slot when provided; takes precedence over `OrgDb`.

- `OrgDb`:

  A Bioconductor `OrgDb` annotation object (e.g., `org.EcK12.eg.db`) for
  GO enrichment. Ignored when `TERM2GENE` is supplied.

- `organism`:

  A KEGG organism code (e.g., `"eco"` for *E. coli* K-12). Requires
  internet access; see
  [`enrichKEGG`](https://rdrr.io/pkg/clusterProfiler/man/enrichKEGG.html).

## GSEA ranking

When `method` includes `"gsea"`, a per-gene score is computed from the
site-level statistics and genes are ranked in decreasing order. Three
metrics are available via `score_metric`:

- `"combined"`:

  (default) \\-\log\_{10}(\text{padj}) \times
  \text{sign}(\Delta\beta)\\; captures both significance and direction.

- `"padj"`:

  \\-\log\_{10}(\text{padj})\\; significance only, direction-agnostic.

- `"delta_beta"`:

  \\\Delta\beta\\ directly; effect-size ranking.

When multiple sites map to the same gene, the score is aggregated across
sites using `gene_score_agg`: `"max"` (default; selects the site with
the largest absolute score, preserving its sign) or `"mean"`.

## See also

[`diffMethyl`](https://carl-stone.github.io/comma/reference/diffMethyl.md),
[`annotateSites`](https://carl-stone.github.io/comma/reference/annotateSites.md),
[`results`](https://carl-stone.github.io/comma/reference/results.md),
[`filterResults`](https://carl-stone.github.io/comma/reference/filterResults.md)

## Examples

``` r
# \donttest{
# Requires clusterProfiler and a custom TERM2GENE mapping
if (requireNamespace("clusterProfiler", quietly = TRUE)) {
  data(comma_example_data)
  dm <- diffMethyl(comma_example_data, formula = ~ condition, mod_type = "6mA")
  ann <- annotateSites(dm, annotation(comma_example_data), type = "overlap")

  # Custom TERM2GENE (works without network access or OrgDb)
  fake_t2g <- data.frame(
    term = c("PATH:01", "PATH:01", "PATH:02"),
    gene = c("geneA",  "geneB",   "geneC")
  )
  res <- enrichMethylation(ann, TERM2GENE = fake_t2g, method = c("ora", "gsea"))
  str(res, max.level = 2)
}
#> 
#> Warning: No significantly differentially methylated genes found (padj <= 0.05 and |delta_beta| >= 0.1).
#> ORA will not be run. Consider relaxing the thresholds.
#> using 'fgsea' for GSEA analysis, please cite Korotkevich et al (2019).
#> preparing geneSet collections...
#> GSEA analysis...
#> no term enriched under specific pvalueCutoff...
#> List of 2
#>  $ go  :List of 2
#>   ..$ ora : NULL
#>   ..$ gsea:Formal class 'gseaResult' [package "DOSE"] with 13 slots
#>  $ kegg:List of 2
#>   ..$ ora : NULL
#>   ..$ gsea: NULL
# }
```
