test_that("parseLocationString handles complex joins", {
  skip_if_not_installed("IRanges")
  parsed <- parseLocationString("complement(join(10..15,20..30))")
  expect_equal(IRanges::start(parsed$ranges), c(10, 20))
  expect_equal(IRanges::end(parsed$ranges), c(15, 30))
  expect_identical(parsed$strand, "-")
})

test_that("parseGenBankSequence removes numbers and whitespace", {
  seq_lines <- c("        1 atgc atgc", "       61 gcta")
  parsed_seq <- parseGenBankSequence(seq_lines)
  expect_equal(parsed_seq, "atgcatgcgcta")
})

test_that("readGenBankFile parses features, qualifiers, and sequence", {
  skip_if_not_installed("GenomicRanges")
  skip_if_not_installed("Biostrings")
  gb_path <- testthat::test_path("fixtures", "mini.gbff")
  parsed <- readGenBankFile(gb_path, include_sequence = TRUE)

  expect_s4_class(parsed$assembly_info, "Seqinfo")
  expect_s4_class(parsed$features, "GRanges")
  expect_equal(seqlevels(parsed$features), "TESTSEQ")

  expect_equal(length(parsed$sequence), 1)
  expect_equal(nchar(as.character(parsed$sequence)), 60)

  gene_feature <- parsed$features[parsed$features$feature_type == "gene"]
  expect_true(all(GenomicRanges::strand(gene_feature) == "-"))
  expect_equal(GenomicRanges::start(gene_feature), 10)
  expect_equal(GenomicRanges::end(gene_feature), 30)
  expect_equal(unique(gene_feature$gene), "geneA")

  cds_features <- parsed$features[parsed$features$feature_type == "CDS"]
  expect_length(cds_features, 2)
  expect_equal(GenomicRanges::start(cds_features), c(10, 20))
  expect_equal(GenomicRanges::end(cds_features), c(15, 30))
  expect_equal(unique(cds_features$product), "Protein A")
  expect_equal(unique(cds_features$note), "example with join")

  split_features <- parsed$features[parsed$features$feature_type == "misc_feature"]
  expect_equal(GenomicRanges::start(split_features), c(60, 70))
  expect_equal(GenomicRanges::end(split_features), c(65, 80))
  expect_equal(unique(split_features$gene), "splitGene")
})
