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

test_that("normalize_annotation_table validates schema and strand values", {
  raw <- data.frame(
    type = "gene",
    id = "geneA",
    seq = "chr1",
    start = 10L,
    end = 20L,
    strand = "+",
    stringsAsFactors = FALSE
  )

  out <- normalize_annotation_table(
    raw,
    feature_type_col = "type",
    feature_id_col = "id",
    seqname_col = "seq",
    start_col = "start",
    end_col = "end",
    strand_col = "strand"
  )

  expect_named(out, c("feature_type", "feature_id", "seqname", "start", "end", "strand"))

  bad <- raw
  bad$strand <- "x"
  expect_error(
    normalize_annotation_table(
      bad,
      feature_type_col = "type",
      feature_id_col = "id",
      seqname_col = "seq",
      start_col = "start",
      end_col = "end",
      strand_col = "strand"
    ),
    "must contain only '\\+', '-', or '\\.'"
  )
})

naive_annotate_sites_with_features <- function(site_table,
                                               feature_table,
                                               site_pos_col = "pos",
                                               site_seqname_col = "seqname",
                                               site_strand_col = NULL,
                                               match_strand = FALSE,
                                               include_unannotated = TRUE) {
  out_rows <- vector("list", length = 0)

  for (i in seq_len(nrow(site_table))) {
    pos <- as.integer(site_table[[site_pos_col]][i])
    seqname <- as.character(site_table[[site_seqname_col]][i])

    hit <- feature_table[
      feature_table$seqname == seqname &
        feature_table$start <= pos &
        feature_table$end >= pos,
      , drop = FALSE
    ]

    if (match_strand) {
      site_strand <- as.character(site_table[[site_strand_col]][i])
      hit <- hit[hit$strand == "." | hit$strand == site_strand, , drop = FALSE]
    }

    if (nrow(hit) == 0) {
      if (isTRUE(include_unannotated)) {
        feature_na <- as.list(setNames(rep(NA_character_, length(names(feature_table))), names(feature_table)))
        out_rows[[length(out_rows) + 1L]] <- c(as.list(site_table[i, , drop = FALSE]), feature_na)
      }
      next
    }

    for (j in seq_len(nrow(hit))) {
      out_rows[[length(out_rows) + 1L]] <- c(as.list(site_table[i, , drop = FALSE]), as.list(hit[j, , drop = FALSE]))
    }
  }

  if (length(out_rows) == 0) {
    out <- site_table[0, , drop = FALSE]
    for (nm in names(feature_table)) {
      out[[nm]] <- character(0)
    }
    return(out)
  }

  out <- do.call(rbind, lapply(out_rows, as.data.frame, stringsAsFactors = FALSE))
  rownames(out) <- NULL
  out
}

test_that("annotate_sites_with_features matches naive output on dense synthetic data", {
  set.seed(1)
  n_contigs <- 20L
  contigs <- paste0("chr", seq_len(n_contigs))

  features <- do.call(
    rbind,
    lapply(contigs, function(chr) {
      starts <- seq(1L, by = 5L, length.out = 80L)
      data.frame(
        feature_type = "window",
        feature_id = paste0(chr, "_", seq_along(starts)),
        seqname = chr,
        start = starts,
        end = starts + 20L,
        strand = rep(c("+", "-", ".", "+"), length.out = length(starts)),
        stringsAsFactors = FALSE
      )
    })
  )

  sites <- data.frame(
    seqname = sample(contigs, 600L, replace = TRUE),
    pos = sample.int(450L, 600L, replace = TRUE),
    strand = sample(c("+", "-"), 600L, replace = TRUE),
    stringsAsFactors = FALSE
  )

  indexed <- annotate_sites_with_features(
    site_table = sites,
    feature_table = features,
    site_strand_col = "strand",
    match_strand = TRUE,
    include_unannotated = TRUE
  )

  naive <- naive_annotate_sites_with_features(
    site_table = sites,
    feature_table = features,
    site_strand_col = "strand",
    match_strand = TRUE,
    include_unannotated = TRUE
  )

  expect_identical(indexed, naive)
})

test_that("annotate_sites_with_features scales better than naive for many contigs", {
  skip_on_cran()
  set.seed(2)

  n_contigs <- 100L
  contigs <- paste0("ctg", seq_len(n_contigs))

  features <- do.call(
    rbind,
    lapply(contigs, function(chr) {
      starts <- seq(1L, by = 4L, length.out = 50L)
      data.frame(
        feature_type = "dense",
        feature_id = paste0(chr, "_", seq_along(starts)),
        seqname = chr,
        start = starts,
        end = starts + 35L,
        strand = ".",
        stringsAsFactors = FALSE
      )
    })
  )

  sites <- data.frame(
    seqname = sample(contigs, 3000L, replace = TRUE),
    pos = sample.int(260L, 3000L, replace = TRUE),
    stringsAsFactors = FALSE
  )

  benchmark_once <- function(expr) {
    system.time(force(expr))[["elapsed"]]
  }

  indexed <- annotate_sites_with_features(
    site_table = sites,
    feature_table = features,
    include_unannotated = FALSE
  )

  naive <- naive_annotate_sites_with_features(
    site_table = sites,
    feature_table = features,
    include_unannotated = FALSE
  )

  time_indexed <- median(replicate(
    3,
    benchmark_once(annotate_sites_with_features(
      site_table = sites,
      feature_table = features,
      include_unannotated = FALSE
    ))
  ))

  time_naive <- median(replicate(
    3,
    benchmark_once(naive_annotate_sites_with_features(
      site_table = sites,
      feature_table = features,
      include_unannotated = FALSE
    ))
  ))

  expect_identical(indexed, naive)
  expect_lte(time_indexed, time_naive)
})
