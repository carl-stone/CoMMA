# commaData: the central data object for the comma package

`commaData` is an S4 class that extends
[`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
to store genome-wide bacterial methylation data from Oxford Nanopore
sequencing. It is the central object accepted and returned by all
`comma` analysis functions.

## Value

An object of class `commaData`. Use
[`commaData`](https://carl-stone.github.io/comma/reference/commaData.md)
to construct instances.

## Details

The class stores methylation data in two assay matrices (accessible via
[`assay`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)):

- `"methylation"`:

  Beta values (proportion of reads called methylated, range 0–1). Sites
  with coverage below the `min_coverage` threshold are stored as `NA`.

- `"coverage"`:

  Integer read depth at each site.

Per-site metadata is in `rowData(object)` and includes at minimum:
`chrom`, `position`, `strand`, `mod_type`, `motif`. The `motif` column
stores the sequence context of each site (e.g., `"GATC"` or `"CCWGG"`)
as extracted from the modkit `mod_code` field. It is `NA` for Dorado and
Megalodon callers, which do not encode per-site motif context.

Per-sample metadata is in `colData(object)` and includes at minimum:
`sample_name`, `condition`, `replicate`.

## Slots

- `genomeInfo`:

  Named integer vector of chromosome sizes c(chromosome name = length in
  bp).

- `annotation`:

  [`GRanges`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
  of genomic features loaded from a GFF3 or BED file. May be an empty
  `GRanges` if no annotation was provided.

- `motifSites`:

  [`GRanges`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
  of all instances of the user-specified sequence motif in the genome
  (e.g., all GATC sites). May be an empty `GRanges` if no motif was
  specified.

## See also

[`commaData`](https://carl-stone.github.io/comma/reference/commaData.md)
for the constructor,
[`methylation`](https://carl-stone.github.io/comma/reference/methylation.md),
[`coverage`](https://rdrr.io/pkg/IRanges/man/coverage-methods.html),
[`sampleInfo`](https://carl-stone.github.io/comma/reference/sampleInfo.md),
[`siteInfo`](https://carl-stone.github.io/comma/reference/siteInfo.md),
[`modTypes`](https://carl-stone.github.io/comma/reference/modTypes.md),
[`annotation`](https://rdrr.io/pkg/BiocGenerics/man/annotation.html) for
accessors.
