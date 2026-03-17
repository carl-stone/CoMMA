# **co**mma: for **m**icrobial **m**ethylation **a**nalysis

The comma package provides a complete toolkit for genome-wide analysis
of bacterial DNA methylation from Oxford Nanopore sequencing data. It
supports the three major modification types (6mA, 5mC, 4mC), handles
input from modkit pileup (primary), Dorado BAM, and Megalodon (legacy)
callers, and provides a unified data container, annotation utilities,
differential methylation testing, and a full set of visualization
functions.

## Value

No return value. This page provides package-level documentation. See
individual function pages for return values.

## Main workflow

1.  **Load data** —
    [`commaData()`](https://carl-stone.github.io/comma/reference/commaData.md)
    constructs the central `commaData` S4 object from per-sample
    methylation files.

2.  **QC** —
    [`methylomeSummary()`](https://carl-stone.github.io/comma/reference/methylomeSummary.md),
    [`coverageDepth()`](https://carl-stone.github.io/comma/reference/coverageDepth.md),
    [`plot_coverage()`](https://carl-stone.github.io/comma/reference/plot_coverage.md),
    [`plot_methylation_distribution()`](https://carl-stone.github.io/comma/reference/plot_methylation_distribution.md),
    and
    [`plot_pca()`](https://carl-stone.github.io/comma/reference/plot_pca.md)
    provide sample-level quality assessment.

3.  **Annotate** —
    [`loadAnnotation()`](https://carl-stone.github.io/comma/reference/loadAnnotation.md)
    imports a GFF3 or BED file;
    [`annotateSites()`](https://carl-stone.github.io/comma/reference/annotateSites.md)
    maps methylation sites to genomic features.

4.  **Visualize** —
    [`plot_genome_track()`](https://carl-stone.github.io/comma/reference/plot_genome_track.md)
    and
    [`plot_metagene()`](https://carl-stone.github.io/comma/reference/plot_metagene.md)
    show methylation in a genomic context.

5.  **Differential methylation** —
    [`diffMethyl()`](https://carl-stone.github.io/comma/reference/diffMethyl.md)
    tests each site;
    [`results()`](https://carl-stone.github.io/comma/reference/results.md)
    and
    [`filterResults()`](https://carl-stone.github.io/comma/reference/filterResults.md)
    extract and filter the results table;
    [`plot_volcano()`](https://carl-stone.github.io/comma/reference/plot_volcano.md)
    and
    [`plot_heatmap()`](https://carl-stone.github.io/comma/reference/plot_heatmap.md)
    visualize the findings.

## Key classes and constructors

- [`commaData`](https://carl-stone.github.io/comma/reference/commaData.md):

  The central S4 data container, extending `SummarizedExperiment`.
  Stores methylation (beta) and coverage matrices, per-site and
  per-sample metadata, genome information, genomic annotation, and motif
  site locations.

## Package options

None. All parameters are passed directly to individual functions.

## References

The modkit pileup format is documented at
<https://nanoporetech.github.io/modkit/>.

## See also

Useful links:

- <https://github.com/carl-stone/CoMMA>

- <https://carl-stone.github.io/CoMMA/>

- Report bugs at <https://github.com/carl-stone/CoMMA/issues>

## Author

**Maintainer**: Carl Stone <carl.j.stone@vanderbilt.edu>
([ORCID](https://orcid.org/0000-0003-0232-5223))
