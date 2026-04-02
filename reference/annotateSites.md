# Annotate methylation sites relative to genomic features

Assigns genomic feature annotations to methylation sites stored in a
[`commaData`](https://carl-stone.github.io/comma/reference/commaData.md)
object using
[`findOverlaps`](https://rdrr.io/pkg/IRanges/man/findOverlaps-methods.html).
Three annotation modes are available: `"overlap"` assigns all
overlapping feature identities to each site, `"proximity"` reports all
features within a distance window and their signed offsets, and
`"metagene"` reports fractional positions within every overlapping
feature.

## Usage

``` r
annotateSites(
  object,
  features = NULL,
  type = c("overlap", "proximity", "metagene"),
  feature_col = "feature_type",
  name_col = "name",
  window = 500L
)
```

## Arguments

- object:

  A
  [`commaData`](https://carl-stone.github.io/comma/reference/commaData.md)
  object.

- features:

  A
  [`GRanges`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
  of genomic features to annotate against. If `NULL` (default), the
  annotation stored in `object` via
  [`annotation`](https://rdrr.io/pkg/BiocGenerics/man/annotation.html)`(object)`
  is used. Must have mcols columns named by `feature_col` and
  `name_col`.

- type:

  Character string specifying the annotation mode. One of:

  `"overlap"`

  :   (default) Each site is assigned all overlapping feature types and
      names. Sites that overlap no feature receive length-0
      `CharacterList` elements.

  `"proximity"`

  :   Each site is assigned all features within `window` bp: their names
      and signed relative positions. Relative positions are 0 when the
      site is inside the feature, negative when upstream (before the
      feature in biological orientation), and positive when downstream.
      Range is \\\[-W, W\]\\ where \\W\\ is `window`. Sites with no
      nearby features receive length-0 elements.

  `"metagene"`

  :   Each site that overlaps a feature is assigned both a fractional
      position (`metagene_frac`, in \\\[0, 1\]\\) and a base-pair offset
      from the biological feature start (`metagene_positions`) for every
      overlapping feature. Strand-aware: for `"-"` strand features,
      position 0 is at the feature end (highest coordinate, the
      biological TSS) and increases toward lower coordinates.
      Non-overlapping sites receive length-0 elements.

- feature_col:

  Character string. Name of the `mcols` column in `features` that
  contains the feature type (e.g., `"feature_type"`). Default:
  `"feature_type"`.

- name_col:

  Character string. Name of the `mcols` column in `features` that
  contains the feature name (e.g., `"name"`). Default: `"name"`.

- window:

  Integer. Window size in base pairs for `type = "proximity"`. All
  features within this distance are returned. Default: `500L`.

## Value

A
[`commaData`](https://carl-stone.github.io/comma/reference/commaData.md)
object identical to `object` except that `rowData` has been extended
with new list-valued annotation columns:

- For `type = "overlap"`::

  `feature_types` (`CharacterList`) and `feature_names`
  (`CharacterList`) — all overlapping feature types and names per site.
  Intergenic sites: `lengths(feature_types) == 0`.

- For `type = "proximity"`::

  `nearby_features` (`CharacterList`) and `rel_positions`
  (`IntegerList`) — all features within `window` bp and their signed
  distances. Sites with none: `lengths(nearby_features) == 0`.

- For `type = "metagene"`::

  `metagene_features` (`CharacterList`), `metagene_frac`
  (`NumericList`), and `metagene_positions` (`IntegerList`) —
  overlapping feature names, fractional positions in \\\[0, 1\]\\, and
  base-pair offsets from the biological feature start. Non-overlapping
  sites: `lengths(metagene_features) == 0`.

## Details

All three modes return every matching feature per site. Results are
stored as
[`CharacterList`](https://rdrr.io/pkg/IRanges/man/AtomicList-class.html),
[`IntegerList`](https://rdrr.io/pkg/IRanges/man/AtomicList-class.html),
or `NumericList` columns in `rowData`. Sites with no overlapping/nearby
features receive length-0 list elements; test for them with
`lengths(col) == 0`.

## Examples

``` r
data(comma_example_data)
# Overlap annotation using built-in annotation
annotated <- annotateSites(comma_example_data)
si <- siteInfo(annotated)
# All overlapping feature types for the first site:
si$feature_types[[1]]
#> [1] "gene"
# Number of sites that overlap at least one feature:
sum(lengths(si$feature_types) > 0)
#> [1] 7
# Intergenic sites:
sum(lengths(si$feature_types) == 0)
#> [1] 293

# Metagene annotation
mg <- annotateSites(comma_example_data, type = "metagene")
si_mg <- siteInfo(mg)
# Fractional and bp positions for the first overlapping site:
si_mg$metagene_frac[[which(lengths(si_mg$metagene_frac) > 0)[1]]]
#> [1] 0.8877756
si_mg$metagene_positions[[which(lengths(si_mg$metagene_positions) > 0)[1]]]
#> [1] 443
```
