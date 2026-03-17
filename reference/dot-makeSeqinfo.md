# Build a Seqinfo object from a named integer vector of chromosome sizes

Build a Seqinfo object from a named integer vector of chromosome sizes

## Usage

``` r
.makeSeqinfo(genome_info, genome_name = NA_character_)
```

## Arguments

- genome_info:

  Named integer vector of chromosome sizes, as returned by
  [`.validateGenomeInfo()`](https://carl-stone.github.io/comma/reference/dot-validateGenomeInfo.md).

- genome_name:

  Optional character string for the genome name.

## Value

A
[`GenomeInfoDb::Seqinfo`](https://rdrr.io/pkg/Seqinfo/man/Seqinfo-class.html)
object, or `NULL` if `genome_info` is `NULL`.
