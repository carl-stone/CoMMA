# Translate and emit collected methylKit warnings as informative summaries

Called after
[`methylKit::calculateDiffMeth()`](https://rdrr.io/pkg/methylKit/man/calculateDiffMeth-methods.html)
to emit one summarized warning per unique warning type rather than
per-site spam. Unknown warning types are re-emitted verbatim so nothing
is silently lost.

## Usage

``` r
.emitMethylKitWarnings(warn_counts)
```

## Arguments

- warn_counts:

  Named list mapping warning message text to integer counts.
