# Vectorized circular genome index

Converts a position that may wrap past the end (or before the start) of
a circular genome into a valid 1-based position.

## Usage

``` r
.circularIndex(position, genome_size)
```

## Arguments

- position:

  Integer vector of genomic positions (1-based).

- genome_size:

  Integer. Total size of the circular genome in base pairs.

## Value

Integer vector of positions in the range `[1, genome_size]`.
