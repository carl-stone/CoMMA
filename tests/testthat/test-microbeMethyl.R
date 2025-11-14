test_that("buildMicrobeMethyl returns MicrobeMethyl object with metadata", {
  bed_data <- data.frame(
    chrom = rep("chrI", 3),
    start = c(10, 20, 30),
    end = c(11, 21, 31),
    name = paste0("site", 1:3),
    score = 0,
    strand = c("+", "+", "-"),
    thickStart = c(10, 20, 30),
    thickEnd = c(11, 21, 31),
    itemRgb = "0",
    coverage = c(50, 60, 70),
    percentMethylation = c(0.9, 0.5, 0.1)
  )
  bed_file <- tempfile(fileext = ".bed")
  write.table(bed_data, file = bed_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

  metadata <- list(condition = "WT", replicate = 1)
  mm <- buildMicrobeMethyl(bed_file, sample_name = "SampleA", metadata = metadata)

  expect_s4_class(mm, "MicrobeMethyl")
  expect_equal(mm@sample_name, "SampleA")
  expect_equal(mm@sample_metadata, metadata)
  expect_equal(mm$coverage, bed_data$coverage)
  expect_equal(mm$percentMethylation, bed_data$percentMethylation)
  expect_s3_class(mm@site_metadata, "data.frame")
  expect_true(all(c("name", "score", "strand", "thickStart", "thickEnd", "itemRgb") %in%
                    colnames(mm@site_metadata)))
})

test_that("MicrobeMethylExperiment stores samples and exposes beta values", {
  mm1 <- MicrobeMethyl(
    assembly = "chrI",
    start = 1:3,
    end = 1:3,
    coverage = c(10, 20, 30),
    percentMethylation = c(0.1, 0.2, 0.3),
    sample_name = "A"
  )
  mm2 <- MicrobeMethyl(
    assembly = "chrI",
    start = 4:6,
    end = 4:6,
    coverage = c(40, 50, 60),
    percentMethylation = c(0.4, 0.5, 0.6),
    sample_name = "B"
  )

  experiment <- MicrobeMethylExperiment(mm1, mm2, experiment_metadata = list(design = "test"))

  expect_s4_class(experiment, "MicrobeMethylExperiment")
  expect_equal(sampleNames(experiment), c("A", "B"))
  beta_values <- getBetaValues(experiment)
  expect_equal(beta_values$A, mm1$percentMethylation)
  expect_equal(beta_values$B, mm2$percentMethylation)
})
