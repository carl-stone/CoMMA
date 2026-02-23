test_that("validate_site_table returns canonical table with beta", {
  input <- data.frame(
    seqname = c("chr1", "chr1"),
    pos = c(10, 11),
    strand = c("+", "-"),
    mod_base = c("6mA", "6mA"),
    n_mod = c(5, 0),
    n_total = c(10, 3),
    sample_id = c("S1", "S1"),
    group = c("WT", "WT"),
    stringsAsFactors = FALSE
  )

  out <- validate_site_table(input)

  expect_true(all(c("seqname", "pos", "strand", "mod_base", "motif", "n_mod", "n_total", "sample_id", "group", "beta") %in% names(out)))
  expect_equal(out$beta, c(0.5, 0))
  expect_equal(unique(out$motif), "GATC")
})

test_that("validate_site_table backfills mod_base and motif defaults for legacy inputs", {
  input <- data.frame(
    seqname = "chr1",
    pos = 10,
    strand = "+",
    n_mod = 5,
    n_total = 10,
    sample_id = "S1",
    group = "WT",
    stringsAsFactors = FALSE
  )

  out <- validate_site_table(input)

  expect_equal(out$mod_base, "6mA")
  expect_equal(out$motif, "GATC")

  out_5mc <- validate_site_table(input, default_mod_base = "5mC", default_motif = "CCWGG")
  expect_equal(out_5mc$mod_base, "5mC")
  expect_equal(out_5mc$motif, "CCWGG")
})

test_that("validate_site_table handles missing counts explicitly", {
  input <- data.frame(
    seqname = "chr1",
    pos = 10,
    strand = "+",
    mod_base = "6mA",
    n_mod = NA_real_,
    n_total = 10,
    sample_id = "S1",
    group = "WT",
    stringsAsFactors = FALSE
  )

  expect_warning(
    out <- validate_site_table(input),
    "setting `beta` to NA"
  )
  expect_true(is.na(out$beta))
})

test_that("validate_site_table errors on missing required columns", {
  bad <- data.frame(
    seqname = "chr1",
    pos = 1,
    strand = "+",
    mod_base = "6mA",
    n_mod = 1,
    n_total = 2,
    sample_id = "S1",
    stringsAsFactors = FALSE
  )

  expect_error(
    validate_site_table(bad),
    "Missing required columns: group"
  )
})

test_that("validate_site_table errors on malformed values", {
  bad <- data.frame(
    seqname = "chr1",
    pos = -1,
    strand = "x",
    mod_base = "6mA",
    n_mod = 2,
    n_total = 1,
    sample_id = "S1",
    group = "WT",
    stringsAsFactors = FALSE
  )

  expect_error(validate_site_table(bad), "`pos` must contain positive integer positions")
})

test_that("convert_modkit_bedmethyl maps and validates modkit output", {
  modkit <- data.frame(
    chrom = c("chr1", "chr1"),
    start = c(0, 9),
    strand = c("+", "-"),
    modified_primary_base = c("6mA", "6mA"),
    n_mod = c(3, 2),
    valid_coverage = c(5, 10),
    stringsAsFactors = FALSE
  )

  out <- convert_modkit_bedmethyl(modkit, sample_id = "S1", group = "WT")

  expect_equal(out$pos, c(1, 10))
  expect_equal(out$n_total, c(5, 10))
  expect_equal(out$sample_id, c("S1", "S1"))
  expect_equal(out$group, c("WT", "WT"))
  expect_equal(out$beta, c(0.6, 0.2))
})

test_that("convert_modkit_bedmethyl errors when required source columns are missing", {
  modkit <- data.frame(
    chrom = "chr1",
    start = 0,
    strand = "+",
    modified_primary_base = "6mA",
    valid_coverage = 5,
    stringsAsFactors = FALSE
  )

  expect_error(
    convert_modkit_bedmethyl(modkit, sample_id = "S1", group = "WT"),
    "missing required columns: n_mod"
  )
})
