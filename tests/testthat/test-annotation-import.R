test_that("read_annotation_gff parses and normalizes canonical feature table", {
  gff_lines <- c(
    "chr1\tsource\tgene\t10\t20\t.\t+\t.\tID=geneA;Name=foo",
    "chr1\tsource\tpromoter\t15\t18\t.\t.\t.\tName=pro1"
  )
  path <- tempfile(fileext = ".gff")
  writeLines(gff_lines, con = path)

  out <- read_annotation_gff(path)

  expect_equal(names(out), c("feature_type", "feature_id", "seqname", "start", "end", "strand", "attributes"))
  expect_equal(out$feature_id[1], "geneA")
  expect_equal(out$feature_id[2], "pro1")
  expect_equal(out$strand, c("+", "."))
})

test_that("read_annotation_bed handles BED coordinates and missing names", {
  bed_lines <- c(
    "chr1\t0\t10\tfeatA\t0\t+",
    "chr1\t10\t15"
  )
  path <- tempfile(fileext = ".bed")
  writeLines(bed_lines, con = path)

  out <- read_annotation_bed(path, feature_type = "misc")

  expect_equal(out$start, c(1L, 11L))
  expect_equal(out$end, c(10L, 15L))
  expect_equal(out$feature_id[1], "featA")
  expect_match(out$feature_id[2], "misc:11-15")
  expect_equal(out$strand, c("+", "."))
})

test_that("annotate_sites_with_features preserves multiple feature hits", {
  sites <- data.frame(
    seqname = c("chr1", "chr1"),
    pos = c(16L, 30L),
    strand = c("+", "+"),
    stringsAsFactors = FALSE
  )

  features <- data.frame(
    feature_type = c("gene", "promoter"),
    feature_id = c("geneA", "proA"),
    seqname = c("chr1", "chr1"),
    start = c(10L, 15L),
    end = c(20L, 18L),
    strand = c("+", "+"),
    stringsAsFactors = FALSE
  )

  out <- annotate_sites_with_features(sites, features)

  expect_equal(sum(out$pos == 16L), 2)
  expect_equal(unique(out$feature_id[out$pos == 16L]), c("geneA", "proA"))
  expect_true(any(out$pos == 30L & is.na(out$feature_id)))
})

test_that("annotate_sites_with_features handles strand matching", {
  sites <- data.frame(
    seqname = "chr1",
    pos = 12L,
    strand = "+",
    stringsAsFactors = FALSE
  )

  features <- data.frame(
    feature_type = c("gene", "gene"),
    feature_id = c("plus_gene", "minus_gene"),
    seqname = c("chr1", "chr1"),
    start = c(10L, 10L),
    end = c(20L, 20L),
    strand = c("+", "-"),
    stringsAsFactors = FALSE
  )

  out <- annotate_sites_with_features(
    sites,
    features,
    site_strand_col = "strand",
    match_strand = TRUE,
    include_unannotated = FALSE
  )

  expect_equal(out$feature_id, "plus_gene")
})
