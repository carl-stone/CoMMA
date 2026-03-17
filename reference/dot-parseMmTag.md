# Parse MM and ML BAM tags into per-base modification calls

Interprets the MM (base modification) and ML (modification likelihood)
tags according to the SAM specification and returns a data frame of read
positions with their modification type and whether they are called
modified.

## Usage

``` r
.parseMmTag(mm_tag, ml_tag, seq_bases)
```

## Arguments

- mm_tag:

  Character string. The MM tag value from the BAM record, e.g.,
  `"A+a?,0,1,3;C+m?,5;"`.

- ml_tag:

  Integer vector. The ML tag values (0–255), parallel to the
  modifications listed in `mm_tag`.

- seq_bases:

  Character string. Read sequence.

## Value

A `data.frame` with columns:

- `read_pos`:

  1-based position in the read.

- `mod_type`:

  Modification type string (e.g., `"6mA"`).

- `is_mod`:

  Logical; `TRUE` when ML probability \> 0.5.

Returns `NULL` on parse failure.
