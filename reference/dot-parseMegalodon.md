# Parse a Megalodon per-read methylation file into a tidy per-site data frame

Reads a Megalodon per-read modification output file, aggregates per-read
calls to per-site beta values and coverage, and returns a tidy data
frame compatible with the
[`commaData`](https://carl-stone.github.io/comma/reference/commaData.md)
constructor. This is an internal function called when
`caller = "megalodon"`.

## Usage

``` r
.parseMegalodon(file, sample_name, mod_type = NULL, min_coverage = 5L)
```

## Arguments

- file:

  Character string. Path to the Megalodon per-read BED file.

- sample_name:

  Character string. Sample name (used in messages).

- mod_type:

  Character string or `NULL`. Modification type to assign to all sites
  (e.g., `"6mA"`). Megalodon files are modification- type-specific, so
  the type cannot be auto-detected from the file alone. If `NULL`,
  defaults to `"6mA"` with a warning.

- min_coverage:

  Integer. Minimum read depth. Default `5`.

## Value

A `data.frame` with columns: `chrom`, `position` (1-based), `strand`,
`mod_type`, `beta`, `coverage`.

## Details

Megalodon per-read output (`modified_bases.5mC.bed` or similar) has the
format:

      chrom  start  end  read_id  score  strand  ...  mod_prob

where `mod_prob` is the per-read probability of modification. This
function aggregates across reads at each site by computing:

- `beta` = mean of per-read probabilities at each site

- `coverage` = number of reads overlapping each site

Sites with `coverage < min_coverage` are dropped.
