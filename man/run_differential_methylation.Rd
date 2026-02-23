% Generated manually because roxygen could not be run in this environment.
\name{run_differential_methylation}
\alias{run_differential_methylation}
\title{Run differential methylation analysis with methylKit}
\usage{
run_differential_methylation(
  site_table,
  mod_base = "6mA",
  motif = "GATC",
  remove_mutated_sites = TRUE,
  mutated_sites = NULL,
  coverage_min = 10,
  qvalue_threshold = 0.05,
  percent_diff_threshold = 10,
  destrand = FALSE,
  context_label = "CpG"
)
}
\arguments{
\item{site_table}{Canonical CoMMA site table with columns at minimum:
`seqname`, `pos`, `strand`, `n_mod`, `n_total`, `sample_id`, and `group`.
`mod_base` and `motif` are supported; missing values are defaulted to
`"6mA"` and `"GATC"` respectively.}

\item{mod_base}{Character scalar specifying which modification type to analyze.
Default `"6mA"` preserves historical behavior. Set to `NULL` to skip
modification-type filtering.}

\item{motif}{Character scalar specifying which motif to analyze. Default
`"GATC"` preserves historical behavior. Set to `NULL` to skip motif filtering.}

\item{remove_mutated_sites}{Logical; whether to remove mutated sites before fitting the model. Default `TRUE`.}

\item{mutated_sites}{Optional mutated-loci specification. Can be `NULL`, numeric positions, or a data frame with `seqname` and `pos` (optional `strand`).}

\item{coverage_min}{Minimum total coverage required to keep a site. Default `10`.}

\item{qvalue_threshold}{Q-value significance threshold. Default `0.05`.}

\item{percent_diff_threshold}{Absolute percent methylation difference threshold. Default `10`.}

\item{destrand}{Logical passed to \code{methylKit::unite()}. Default `FALSE`.}

\item{context_label}{Character scalar metadata label stored in methylKit objects via\n\code{methylRaw(..., context = context_label)}. This should match the analyzed\nmodification/motif context in bacterial workflows. Default `"CpG"` to\npreserve historical behavior.}
}
\value{
A named list with:
\itemize{
\item `result_table`: tidy differential methylation table containing `seqname`, `pos`, `strand`, `pvalue`, `qvalue`, `percent_diff`, `significant`.
\item `plot_data`: list with `volcano_data` and `beta_by_sample`.
\item `methylkit_objects`: list with methylKit intermediate objects (`methyl_raw`, `normalized`, `united`, `diff`).
}
}
\description{
High-level differential methylation wrapper that applies CoMMA defaults:
coverage filtering, median normalization, methylKit logistic regression, and
SLIM correction.
}
