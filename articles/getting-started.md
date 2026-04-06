# Getting Started with comma

## Introduction

`comma` (**Co**mparative **M**ethylomics for **M**icrobial **A**nalysis)
is an R package for genome-wide analysis of bacterial DNA methylation
from Oxford Nanopore sequencing data. It supports three modification
types — N6-methyladenine (6mA), 5-methylcytosine (5mC), and
N4-methylcytosine (4mC) — in a single, unified data container. This
vignette walks through the complete analysis workflow using the built-in
`comma_example_data` synthetic dataset.

The typical `comma` workflow has five steps:

1.  **Load** per-sample methylation files into a `commaData` object.
2.  **QC** the data (coverage, beta distributions, PCA).
3.  **Annotate** sites relative to genomic features.
4.  **Test** for differential methylation between conditions.
5.  **Visualize** the results.

## Installation

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("comma")
```

``` r
library(comma)
```

## The `commaData` Object

`commaData` extends `SummarizedExperiment` and is the central data
container in `comma`. It stores:

- **methylation** — a sites × samples matrix of beta values (0–1).
- **coverage** — a sites × samples matrix of read depths.
- **rowData** — per-site metadata: chromosome, position, strand,
  mod_type.
- **colData** — per-sample metadata: sample_name, condition, replicate.
- **genomeInfo** — chromosome names and lengths.
- **annotation** — genomic features as a `GRanges` object.
- **motifSites** — motif instances as a `GRanges` object.

The built-in `comma_example_data` contains 300 synthetic methylation
sites (200 × 6mA, 100 × 5mC) on a simulated 100 kb chromosome across
three samples: two controls (`ctrl_1`, `ctrl_2`) and one treatment
(`treat_1`).

``` r
data(comma_example_data)
comma_example_data
#> class: commaData
#> sites: 588 | samples: 6 
#> mod types: 5mC, 6mA 
#> motifs: CCWGG, GATC 
#> mod contexts: 5mC_CCWGG, 6mA_GATC 
#> conditions: control, treatment 
#> genome: 1 chromosome (100,000 bp total) 
#> annotation: 5 features 
#> motif sites: none
```

``` r
# Modification types present
modTypes(comma_example_data)
#> [1] "5mC" "6mA"

# Per-sample metadata
sampleInfo(comma_example_data)
#>         sample_name condition replicate caller
#> ctrl_1       ctrl_1   control         1 modkit
#> ctrl_2       ctrl_2   control         2 modkit
#> ctrl_3       ctrl_3   control         3 modkit
#> treat_1     treat_1 treatment         1 modkit
#> treat_2     treat_2 treatment         2 modkit
#> treat_3     treat_3 treatment         3 modkit

# Matrix dimensions: sites × samples
dim(methylation(comma_example_data))
#> [1] 588   6
```

## Exploring the Methylome

### Summary Statistics

[`methylomeSummary()`](https://carl-stone.github.io/comma/reference/methylomeSummary.md)
returns a tidy data frame with per-sample distribution statistics:

``` r
ms <- methylomeSummary(comma_example_data)
ms[, c("sample_name", "condition", "mean_beta", "median_beta", "n_covered")]
#>   sample_name condition mean_beta median_beta n_covered
#> 1      ctrl_1   control 0.8654839   0.8929141       588
#> 2      ctrl_2   control 0.8705692   0.8959600       588
#> 3      ctrl_3   control 0.8638033   0.8918851       588
#> 4     treat_1 treatment 0.8357998   0.8864176       588
#> 5     treat_2 treatment 0.8369054   0.8893089       588
#> 6     treat_3 treatment 0.8388398   0.8866568       588
```

### Coverage QC

[`plot_coverage()`](https://carl-stone.github.io/comma/reference/plot_coverage.md)
shows the distribution of sequencing depth per site, per sample.
Consistent coverage across samples is an important quality indicator.

``` r
plot_coverage(comma_example_data)
```

![Coverage depth distribution per
sample.](getting-started_files/figure-html/plot-coverage-1.png)

Coverage depth distribution per sample.

### Beta Value Distributions

[`plot_methylation_distribution()`](https://carl-stone.github.io/comma/reference/plot_methylation_distribution.md)
plots the density of methylation levels for each sample. Bacterial
genomes often show a bimodal distribution (sites are either fully
methylated or unmethylated).

``` r
plot_methylation_distribution(comma_example_data)
```

![Methylation beta value density per sample, faceted by modification
type.](getting-started_files/figure-html/plot-dist-1.png)

Methylation beta value density per sample, faceted by modification type.

Restrict to a single modification type:

``` r
plot_methylation_distribution(comma_example_data, mod_type = "6mA")
```

![Beta value density for 6mA sites
only.](getting-started_files/figure-html/plot-dist-6ma-1.png)

Beta value density for 6mA sites only.

### PCA for Sample-Level QC

[`plot_pca()`](https://carl-stone.github.io/comma/reference/plot_pca.md)
performs PCA on per-sample methylation profiles. Samples from the same
condition should cluster together. Internally, beta values are converted
to M-values via
[`mValues()`](https://carl-stone.github.io/comma/reference/mValues.md)
before PCA, which stabilizes variance across sites near 0 or 1.

``` r
plot_pca(comma_example_data, color_by = "condition")
```

![PCA of methylation profiles colored by
condition.](getting-started_files/figure-html/plot-pca-1.png)

PCA of methylation profiles colored by condition.

To retrieve the underlying scores for custom plotting, use
`return_data = TRUE`. The result is a `data.frame` with `PC1`, `PC2`,
and all sample metadata columns; the percentage of variance explained by
each PC is stored in `attr(result, "percentVar")`.

``` r
pca_df <- plot_pca(comma_example_data, return_data = TRUE)
attr(pca_df, "percentVar")  # variance explained by PC1 and PC2
#>  PC1  PC2 
#> 36.3 17.4
```

## Annotating Sites

[`annotateSites()`](https://carl-stone.github.io/comma/reference/annotateSites.md)
maps methylation sites to genomic features, always computing four
parallel list columns in `rowData`:

- `feature_types` — GFF3 feature type for each association (e.g.,
  `"gene"`, `"promoter"`).
- `feature_names` — feature name for each association.
- `rel_position` — signed distance from the feature (0 = inside;
  negative = upstream; positive = downstream, strand-aware).
- `frac_position` — normalized position within the feature (\[0, 1\];
  `NA` for sites outside).

Use the `keep` argument to filter which associations are retained:
`"all"` (default), `"overlap"` (only overlapping features),
`"proximity"` (retains `rel_position`, drops `frac_position`), or
`"metagene"` (only overlapping features with `frac_position`).

``` r
annotated <- annotateSites(comma_example_data)
si <- siteInfo(annotated)

# Proportion of sites overlapping at least one annotated feature
mean(lengths(si$feature_names) > 0)
#> [1] 0.03401361
```

[`plot_metagene()`](https://carl-stone.github.io/comma/reference/plot_metagene.md)
visualizes the average methylation profile across gene bodies:

``` r
plot_metagene(comma_example_data, feature = "gene")
```

![Mean methylation profile across gene bodies (TSS to
TTS).](getting-started_files/figure-html/plot-metagene-1.png)

Mean methylation profile across gene bodies (TSS to TTS).

## Genome Track Visualization

[`plot_genome_track()`](https://carl-stone.github.io/comma/reference/plot_genome_track.md)
produces a genome browser–style plot of methylation along a chromosome
region:

``` r
plot_genome_track(comma_example_data, chromosome = "chr_sim",
                  start = 1L, end = 50000L, mod_type = "6mA")
```

![Genome track for the first 50 kb of
chr_sim.](getting-started_files/figure-html/plot-track-1.png)

Genome track for the first 50 kb of chr_sim.

## Differential Methylation

[`diffMethyl()`](https://carl-stone.github.io/comma/reference/diffMethyl.md)
tests each site for differential methylation between conditions. It is
modeled on DESeq2’s workflow: pass a `commaData` object and a design
formula, and receive back the same object with statistical results in
`rowData`.

``` r
cd_dm <- diffMethyl(comma_example_data, formula = ~ condition,
                    mod_type = "6mA")
cd_dm
#> class: commaData
#> sites: 588 | samples: 6 
#> mod types: 5mC, 6mA 
#> motifs: CCWGG, GATC 
#> mod contexts: 5mC_CCWGG, 6mA_GATC 
#> conditions: control, treatment 
#> genome: 1 chromosome (100,000 bp total) 
#> annotation: 5 features 
#> motif sites: none
```

Extract the results as a tidy data frame:

``` r
res <- results(cd_dm)
# Top sites by adjusted p-value
head(res[order(res$dm_padj),
         c("chrom", "position", "mod_type", "dm_delta_beta", "dm_padj")])
#>                            chrom position mod_type dm_delta_beta      dm_padj
#> chr_sim:50176:-:6mA:GATC chr_sim    50176      6mA    -0.7336497 1.849154e-75
#> chr_sim:70003:-:6mA:GATC chr_sim    70003      6mA    -0.7050844 3.896483e-68
#> chr_sim:63550:+:6mA:GATC chr_sim    63550      6mA    -0.7799241 5.006897e-66
#> chr_sim:61440:+:6mA:GATC chr_sim    61440      6mA    -0.7090099 1.178364e-64
#> chr_sim:86016:+:6mA:GATC chr_sim    86016      6mA    -0.6743832 3.661541e-62
#> chr_sim:2180:-:6mA:GATC  chr_sim     2180      6mA    -0.7543758 4.024452e-60
```

Filter to significant sites (padj \< 0.05, \|Δβ\| ≥ 0.2):

``` r
sig <- filterResults(cd_dm, padj = 0.05, delta_beta = 0.2)
cat("Significant sites:", nrow(sig), "\n")
#> Significant sites: 31
```

### Volcano Plot

[`plot_volcano()`](https://carl-stone.github.io/comma/reference/plot_volcano.md)
displays the differential methylation landscape. Sites are colored by
direction and significance:

``` r
plot_volcano(res)
```

![Volcano plot: effect size (Δβ) vs. significance (–log10
padj).](getting-started_files/figure-html/plot-volcano-1.png)

Volcano plot: effect size (Δβ) vs. significance (–log10 padj).

### Heatmap of Top Sites

[`plot_heatmap()`](https://carl-stone.github.io/comma/reference/plot_heatmap.md)
shows methylation beta values for the top differentially methylated
sites:

``` r
plot_heatmap(res, cd_dm, n_sites = 30L)
```

![Heatmap of top 30 differentially methylated 6mA
sites.](getting-started_files/figure-html/plot-heatmap-1.png)

Heatmap of top 30 differentially methylated 6mA sites.

## Session Information

``` r
sessionInfo()
#> R version 4.5.3 (2026-03-11)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.4 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices datasets  utils     methods   base     
#> 
#> other attached packages:
#> [1] comma_0.7.2.9000 BiocStyle_2.38.0
#> 
#> loaded via a namespace (and not attached):
#>   [1] bitops_1.0-9                rlang_1.1.7                
#>   [3] magrittr_2.0.5              matrixStats_1.5.0          
#>   [5] compiler_4.5.3              mgcv_1.9-4                 
#>   [7] systemfonts_1.3.2           vctrs_0.7.2                
#>   [9] reshape2_1.4.5              stringr_1.6.0              
#>  [11] pkgconfig_2.0.3             crayon_1.5.3               
#>  [13] fastmap_1.2.0               XVector_0.50.0             
#>  [15] labeling_0.4.3              Rsamtools_2.26.0           
#>  [17] rmarkdown_2.31              UCSC.utils_1.6.1           
#>  [19] ragg_1.5.2                  xfun_0.57                  
#>  [21] cachem_1.1.0                cigarillo_1.0.0            
#>  [23] GenomeInfoDb_1.46.2         jsonlite_2.0.0             
#>  [25] DelayedArray_0.36.1         BiocParallel_1.44.0        
#>  [27] parallel_4.5.3              R6_2.6.1                   
#>  [29] bslib_0.10.0                stringi_1.8.7              
#>  [31] RColorBrewer_1.1-3          limma_3.66.0               
#>  [33] rtracklayer_1.70.1          GenomicRanges_1.62.1       
#>  [35] jquerylib_0.1.4             numDeriv_2016.8-1.1        
#>  [37] Rcpp_1.1.1                  Seqinfo_1.0.0              
#>  [39] bookdown_0.46               SummarizedExperiment_1.40.0
#>  [41] knitr_1.51                  zoo_1.8-15                 
#>  [43] R.utils_2.13.0              IRanges_2.44.0             
#>  [45] Matrix_1.7-4                splines_4.5.3              
#>  [47] tidyselect_1.2.1            qvalue_2.42.0              
#>  [49] abind_1.4-8                 yaml_2.3.12                
#>  [51] codetools_0.2-20            curl_7.0.0                 
#>  [53] lattice_0.22-9              tibble_3.3.1               
#>  [55] plyr_1.8.9                  Biobase_2.70.0             
#>  [57] withr_3.0.2                 S7_0.2.1                   
#>  [59] coda_0.19-4.1               evaluate_1.0.5             
#>  [61] desc_1.4.3                  mclust_6.1.2               
#>  [63] Biostrings_2.78.0           pillar_1.11.1              
#>  [65] BiocManager_1.30.27         MatrixGenerics_1.22.0      
#>  [67] renv_1.1.8                  stats4_4.5.3               
#>  [69] generics_0.1.4              RCurl_1.98-1.18            
#>  [71] emdbook_1.3.14              S4Vectors_0.48.1           
#>  [73] ggplot2_4.0.2               scales_1.4.0               
#>  [75] gtools_3.9.5                glue_1.8.0                 
#>  [77] tools_4.5.3                 BiocIO_1.20.0              
#>  [79] data.table_1.18.2.1         GenomicAlignments_1.46.0   
#>  [81] fs_2.0.1                    mvtnorm_1.3-6              
#>  [83] XML_3.99-0.23               grid_4.5.3                 
#>  [85] bbmle_1.0.25.1              bdsmatrix_1.3-7            
#>  [87] nlme_3.1-168                patchwork_1.3.2            
#>  [89] restfulr_0.0.16             cli_3.6.5                  
#>  [91] textshaping_1.0.5           fastseg_1.56.0             
#>  [93] S4Arrays_1.10.1             methylKit_1.36.0           
#>  [95] dplyr_1.2.1                 gtable_0.3.6               
#>  [97] R.methodsS3_1.8.2           sass_0.4.10                
#>  [99] digest_0.6.39               BiocGenerics_0.56.0        
#> [101] SparseArray_1.10.10         rjson_0.2.23               
#> [103] htmlwidgets_1.6.4           farver_2.1.2               
#> [105] htmltools_0.5.9             pkgdown_2.2.0              
#> [107] R.oo_1.27.1                 lifecycle_1.0.5            
#> [109] httr_1.4.8                  statmod_1.5.1              
#> [111] MASS_7.3-65
```
