# Package index

## All functions

- [`annotateSites()`](https://carl-stone.github.io/comma/reference/annotateSites.md)
  : Annotate methylation sites relative to genomic features

- [`annotation(`*`<commaData>`*`)`](https://carl-stone.github.io/comma/reference/annotation-commaData-method.md)
  : Accessor for genomic feature annotation

- [`comma-package`](https://carl-stone.github.io/comma/reference/comma-package.md)
  [`comma`](https://carl-stone.github.io/comma/reference/comma-package.md)
  :

  **co**mma: for **m**icrobial **m**ethylation **a**nalysis

- [`commaData-class`](https://carl-stone.github.io/comma/reference/commaData-class.md)
  : commaData: the central data object for the comma package

- [`commaData()`](https://carl-stone.github.io/comma/reference/commaData.md)
  : Create a commaData object from methylation calling output files

- [`comma_example_data`](https://carl-stone.github.io/comma/reference/comma_example_data.md)
  : Synthetic example methylation dataset for the comma package

- [`coverage(`*`<commaData>`*`)`](https://carl-stone.github.io/comma/reference/coverage-commaData-method.md)
  : Accessor for the sequencing coverage (read depth) matrix

- [`coverageDepth()`](https://carl-stone.github.io/comma/reference/coverageDepth.md)
  : Windowed sequencing depth across the genome

- [`diffMethyl()`](https://carl-stone.github.io/comma/reference/diffMethyl.md)
  : Identify differentially methylated sites between conditions

- [`enrichMethylation()`](https://carl-stone.github.io/comma/reference/enrichMethylation.md)
  : Gene set enrichment analysis of differential methylation results

- [`filterResults()`](https://carl-stone.github.io/comma/reference/filterResults.md)
  : Filter differential methylation results by significance thresholds

- [`findMotifSites()`](https://carl-stone.github.io/comma/reference/findMotifSites.md)
  : Find all instances of a sequence motif in a genome

- [`genome(`*`<commaData>`*`)`](https://carl-stone.github.io/comma/reference/genome-commaData-method.md)
  : Accessor for genome size information

- [`loadAnnotation()`](https://carl-stone.github.io/comma/reference/loadAnnotation.md)
  : Load genomic feature annotations from a GFF3 or BED file

- [`mValues()`](https://carl-stone.github.io/comma/reference/mValues.md)
  : Compute M-values from a commaData object

- [`methylation()`](https://carl-stone.github.io/comma/reference/methylation.md)
  : Accessor for the methylation (beta value) matrix

- [`methylomeSummary()`](https://carl-stone.github.io/comma/reference/methylomeSummary.md)
  : Summarize per-sample methylation and coverage distributions

- [`modContexts()`](https://carl-stone.github.io/comma/reference/modContexts.md)
  : Return the modification contexts present in a commaData object

- [`modTypes()`](https://carl-stone.github.io/comma/reference/modTypes.md)
  : Return the modification types present in a commaData object

- [`motifSites()`](https://carl-stone.github.io/comma/reference/motifSites.md)
  : Accessor for motif site positions

- [`motifs()`](https://carl-stone.github.io/comma/reference/motifs.md) :
  Accessor for sequence context motifs present in a commaData object

- [`plot_coverage()`](https://carl-stone.github.io/comma/reference/plot_coverage.md)
  : Plot coverage depth distribution

- [`plot_genome_track()`](https://carl-stone.github.io/comma/reference/plot_genome_track.md)
  : Genome browser-style methylation track plot

- [`plot_heatmap()`](https://carl-stone.github.io/comma/reference/plot_heatmap.md)
  : Heatmap of top differentially methylated sites

- [`plot_metagene()`](https://carl-stone.github.io/comma/reference/plot_metagene.md)
  : Metagene plot of methylation across genomic features

- [`plot_methylation_distribution()`](https://carl-stone.github.io/comma/reference/plot_methylation_distribution.md)
  : Plot methylation beta value distributions

- [`plot_pca()`](https://carl-stone.github.io/comma/reference/plot_pca.md)
  : PCA of methylation profiles

- [`plot_tss_profile()`](https://carl-stone.github.io/comma/reference/plot_tss_profile.md)
  : TSS-centered methylation profile

- [`plot_volcano()`](https://carl-stone.github.io/comma/reference/plot_volcano.md)
  : Volcano plot for differential methylation results

- [`results()`](https://carl-stone.github.io/comma/reference/results.md)
  : Extract differential methylation results as a tidy data frame

- [`sampleInfo()`](https://carl-stone.github.io/comma/reference/sampleInfo.md)
  : Accessor for per-sample metadata

- [`siteInfo()`](https://carl-stone.github.io/comma/reference/siteInfo.md)
  : Accessor for per-site metadata

- [`slidingWindow()`](https://carl-stone.github.io/comma/reference/slidingWindow.md)
  : Sliding window methylation summary along the genome

- [`` `[`( ``*`<commaData>`*`,`*`<ANY>`*`,`*`<ANY>`*`,`*`<ANY>`*`)`](https://carl-stone.github.io/comma/reference/sub-commaData-ANY-ANY-ANY-method.md)
  : Subset a commaData object by sites and/or samples

- [`subset()`](https://carl-stone.github.io/comma/reference/subset.md) :
  Subset a commaData object by condition, modification type, or
  chromosome

- [`varianceByDepth()`](https://carl-stone.github.io/comma/reference/varianceByDepth.md)
  : Methylation variance as a function of sequencing depth

- [`writeBED()`](https://carl-stone.github.io/comma/reference/writeBED.md)
  : Export methylation data as a BED file
