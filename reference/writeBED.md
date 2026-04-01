# Export methylation data as a BED file

Writes per-site methylation beta values for a single sample from a
[`commaData`](https://carl-stone.github.io/comma/reference/commaData.md)
object to a 9-column BED file (BED9 format), suitable for visualisation
in IGV, UCSC Genome Browser, or other genome browsers that support the
`itemRGB` field.

## Usage

``` r
writeBED(
  object,
  file,
  sample,
  mod_type = NULL,
  rgb_scale = TRUE,
  track_name = sample,
  track_description = "methylation beta values"
)
```

## Arguments

- object:

  A
  [`commaData`](https://carl-stone.github.io/comma/reference/commaData.md)
  object.

- file:

  Character string. Path to the output BED file. The file will be
  created or overwritten.

- sample:

  Character string. Name of the sample to export. Must match a column
  name in `methylation(object)`.

- mod_type:

  Character string or `NULL`. If provided, only sites of the specified
  modification type are written (e.g., `"6mA"`). If `NULL` (default),
  all sites are written.

- rgb_scale:

  Logical. If `TRUE` (default), an `itemRGB` column is added based on
  the methylation score using a blue-to-red gradient (low = blue, high =
  red). If `FALSE`, `itemRGB` is set to `"0,0,0"` for all sites.

- track_name:

  Character string. Name shown in the genome browser track header.
  Defaults to `sample`.

- track_description:

  Character string. Description shown in the track header. Defaults to
  `"methylation beta values"`.

## Value

Invisibly returns the path to the written file (`file`).

## Details

The BED score field (column 5) contains the beta value multiplied by
1000 and rounded to the nearest integer (range 0–1000), which is the
standard convention for methylation BED files. The `itemRGB` colour
gradient transitions as follows: score ≤ 200 = blue (0,0,255), score ≤
400 = blue-purple (83,0,172), score ≤ 600 = purple (167,0,85), score ≤
800 = red-purple (222,0,28), score ≤ 1000 = red (250,0,0).

Sites with `NA` methylation (below the coverage threshold) are excluded
from the output.

## See also

[`methylation`](https://carl-stone.github.io/comma/reference/methylation.md),
[`siteInfo`](https://carl-stone.github.io/comma/reference/siteInfo.md)

## Examples

``` r
data(comma_example_data)
tmp <- tempfile(fileext = ".bed")
writeBED(comma_example_data, file = tmp, sample = "ctrl_1", mod_type = "6mA")
#> Warning: missing package slot (comma) in object of class “commaData” (package info added)

# \donttest{
# Write to a permanent file
writeBED(comma_example_data,
         file    = "ctrl_1_methylation.bed",
         sample  = "ctrl_1",
         mod_type = "6mA")
#> Warning: missing package slot (comma) in object of class “commaData” (package info added)
# }
```
