make_site_table <- function() {
  data.frame(
    seqname = rep("chr1", 8),
    pos = c(100, 101, 100, 101, 100, 101, 100, 101),
    strand = rep("+", 8),
    mod_base = rep("6mA", 8),
    motif = rep("GATC", 8),
    n_mod = c(5, 1, 6, 2, 9, 4, 8, 3),
    n_total = c(10, 9, 10, 12, 14, 11, 13, 10),
    sample_id = rep(c("S1", "S2", "S3", "S4"), each = 2),
    group = rep(c("A", "A", "B", "B"), each = 2),
    mutated = c(FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE),
    stringsAsFactors = FALSE
  )
}

test_that("run_differential_methylation enforces default filtering, normalization, model and thresholds", {
  input <- make_site_table()

  seen <- new.env(parent = emptyenv())

  local_mocked_bindings(
    .site_table_to_methylkit = function(site_table) {
      seen$filtered <- site_table
      "raw_obj"
    },
    .normalize_methylkit_coverage = function(mk_raw) {
      seen$normalize_input <- mk_raw
      "norm_obj"
    },
    unite = function(mk_norm, destrand = FALSE) {
      seen$unite_input <- mk_norm
      seen$destrand <- destrand
      "united_obj"
    },
    .fit_methylkit_diff = function(mk_united) {
      seen$fit_input <- mk_united
      "diff_obj"
    },
    .methylkit_diff_to_df = function(mk_diff) {
      seen$diff_input <- mk_diff
      data.frame(
        chr = c("chr1", "chr1"),
        start = c(100, 101),
        strand = c("+", "+"),
        pvalue = c(0.001, 0.2),
        qvalue = c(0.04, 0.04),
        meth.diff = c(12, 8),
        stringsAsFactors = FALSE
      )
    }
  )

  out <- run_differential_methylation(input)

  # defaults preserve 6mA/GATC selection
  expect_true(all(seen$filtered$mod_base == "6mA"))
  expect_true(all(seen$filtered$motif == "GATC"))

  # coverage_min default is 10 and mutated sites are removed by default
  expect_true(all(seen$filtered$n_total >= 10))
  expect_false(any(seen$filtered$mutated))

  # defaults pass through to the methylKit pipeline
  expect_identical(seen$normalize_input, "raw_obj")
  expect_identical(seen$unite_input, "norm_obj")
  expect_identical(seen$fit_input, "united_obj")
  expect_identical(seen$diff_input, "diff_obj")
  expect_false(seen$destrand)

  # default thresholds: q < 0.05 and abs(percent_diff) > 10
  expect_equal(out$result_table$significant, c(TRUE, FALSE))

  # stable output schema
  expect_named(
    out$result_table,
    c("seqname", "pos", "strand", "pvalue", "qvalue", "percent_diff", "significant")
  )
  expect_true(all(c("volcano_data", "beta_by_sample") %in% names(out$plot_data)))
  expect_named(
    out$plot_data$beta_by_sample,
    c("seqname", "pos", "strand", "mod_base", "sample_id", "group", "beta")
  )
})

test_that("run_differential_methylation can keep mutated sites when requested", {
  input <- make_site_table()
  seen <- new.env(parent = emptyenv())

  local_mocked_bindings(
    .site_table_to_methylkit = function(site_table) {
      seen$filtered <- site_table
      "raw_obj"
    },
    .normalize_methylkit_coverage = function(mk_raw) "norm_obj",
    unite = function(mk_norm, destrand = FALSE) "united_obj",
    .fit_methylkit_diff = function(mk_united) "diff_obj",
    .methylkit_diff_to_df = function(mk_diff) {
      data.frame(
        chr = "chr1",
        start = 100,
        strand = "+",
        pvalue = 0.1,
        qvalue = 0.1,
        meth.diff = 2,
        stringsAsFactors = FALSE
      )
    }
  )

  run_differential_methylation(input, remove_mutated_sites = FALSE)
  expect_true(any(seen$filtered$mutated))
})

test_that("run_differential_methylation supports explicit mutated_sites and custom thresholds", {
  input <- make_site_table()

  local_mocked_bindings(
    .site_table_to_methylkit = function(site_table) {
      expect_false(any(site_table$pos %in% 100))
      "raw_obj"
    },
    .normalize_methylkit_coverage = function(mk_raw) "norm_obj",
    unite = function(mk_norm, destrand = FALSE) "united_obj",
    .fit_methylkit_diff = function(mk_united) "diff_obj",
    .methylkit_diff_to_df = function(mk_diff) {
      data.frame(
        chr = "chr1",
        start = 101,
        strand = "+",
        pvalue = 0.01,
        qvalue = 0.049,
        meth.diff = 9,
        stringsAsFactors = FALSE
      )
    }
  )

  out <- run_differential_methylation(
    input,
    mutated_sites = 100,
    qvalue_threshold = 0.05,
    percent_diff_threshold = 8
  )
  expect_true(out$result_table$significant)
})

test_that("run_differential_methylation supports non-6mA data and optional motif filtering", {
  input <- make_site_table()
  input$mod_base <- "5mC"
  input$motif <- "CCWGG"

  seen <- new.env(parent = emptyenv())

  local_mocked_bindings(
    .site_table_to_methylkit = function(site_table) {
      seen$filtered <- site_table
      "raw_obj"
    },
    .normalize_methylkit_coverage = function(mk_raw) "norm_obj",
    unite = function(mk_norm, destrand = FALSE) "united_obj",
    .fit_methylkit_diff = function(mk_united) "diff_obj",
    .methylkit_diff_to_df = function(mk_diff) {
      data.frame(
        chr = "chr1",
        start = 100,
        strand = "+",
        pvalue = 0.02,
        qvalue = 0.02,
        meth.diff = 20,
        stringsAsFactors = FALSE
      )
    }
  )

  run_differential_methylation(input, mod_base = "5mC", motif = "CCWGG")
  expect_true(all(seen$filtered$mod_base == "5mC"))
  expect_true(all(seen$filtered$motif == "CCWGG"))

  run_differential_methylation(input, mod_base = "5mC", motif = NULL)
  expect_true(all(seen$filtered$mod_base == "5mC"))
})
