# Map read positions to reference (genomic) positions via CIGAR string

Returns a vector of length equal to the read sequence, where element `i`
is the 1-based reference position corresponding to read base `i`, or
`NA` for soft-clipped or inserted bases.

## Usage

``` r
.cigarToRefPos(cigar_str, ref_start, seq_bases)
```

## Arguments

- cigar_str:

  Character string. CIGAR string (e.g., `"5M2I3M"`).

- ref_start:

  Integer. 1-based leftmost reference position of the alignment.

- seq_bases:

  Character string of read sequence bases.

## Value

Integer vector of length `nchar(seq_bases)`, or `NULL` on parse failure.
