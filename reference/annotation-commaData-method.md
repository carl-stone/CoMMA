# Accessor for genomic feature annotation

Returns the
[`GRanges`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html) of
genomic features stored in a
[`commaData`](https://carl-stone.github.io/comma/reference/commaData.md)
object. This is the annotation loaded from a GFF3 or BED file at
construction time.

## Usage

``` r
# S4 method for class 'commaData'
annotation(object)
```

## Arguments

- object:

  A `commaData` object.

## Value

A [`GRanges`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
object. May be empty (length 0) if no annotation was provided when
creating the object.

## Examples

``` r
data(comma_example_data)
annotation(comma_example_data)
#> GRanges object with 5 ranges and 2 metadata columns:
#>       seqnames    ranges strand | feature_type        name
#>          <Rle> <IRanges>  <Rle> |  <character> <character>
#>   [1]  chr_sim     1-500      + |         gene       geneA
#>   [2]  chr_sim  600-1200      + |         gene       geneB
#>   [3]  chr_sim 1400-2000      - |         gene       geneC
#>   [4]  chr_sim 2500-3500      + |         rRNA       geneD
#>   [5]  chr_sim 4000-5000      - |         tRNA       geneE
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```
