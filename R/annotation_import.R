#' Normalize parsed annotations to CoMMA's canonical feature table
#'
#' CoMMA uses a canonical feature table for genomic annotations with the
#' following required columns: `feature_type`, `feature_id`, `seqname`, `start`,
#' `end`, and `strand`. An optional `attributes` column can store parser-specific
#' metadata (for example GFF attributes or extra BED columns).
#'
#' This keeps core annotation import organism-agnostic. Organism-specific
#' resources (for example RegulonDB adapters) should be implemented as optional
#' add-ons that emit this canonical schema.
#'
#' @param annotation_df A parsed annotation `data.frame`.
#' @param feature_type_col Column name containing feature types.
#' @param feature_id_col Column name containing feature identifiers.
#' @param seqname_col Column name containing sequence/chromosome names.
#' @param start_col Column name containing 1-based feature start positions.
#' @param end_col Column name containing 1-based inclusive feature end positions.
#' @param strand_col Column name containing strand (`"+"`, `"-"`, or `"."`).
#' @param attributes_col Optional column containing free-form attributes.
#'
#' @return A normalized canonical feature table.
#' @export
normalize_annotation_table <- function(annotation_df,
                                       feature_type_col,
                                       feature_id_col,
                                       seqname_col,
                                       start_col,
                                       end_col,
                                       strand_col = NULL,
                                       attributes_col = NULL) {
  if (!is.data.frame(annotation_df)) {
    stop("`annotation_df` must be a data.frame.", call. = FALSE)
  }

  required <- c(feature_type_col, feature_id_col, seqname_col, start_col, end_col)
  missing_cols <- setdiff(required, names(annotation_df))
  if (length(missing_cols) > 0) {
    stop(
      paste0("Missing required annotation columns: ", paste(missing_cols, collapse = ", "), "."),
      call. = FALSE
    )
  }

  out <- data.frame(
    feature_type = as.character(annotation_df[[feature_type_col]]),
    feature_id = as.character(annotation_df[[feature_id_col]]),
    seqname = as.character(annotation_df[[seqname_col]]),
    start = as.integer(annotation_df[[start_col]]),
    end = as.integer(annotation_df[[end_col]]),
    stringsAsFactors = FALSE
  )

  if (!is.null(strand_col)) {
    if (!(strand_col %in% names(annotation_df))) {
      stop("`strand_col` is not present in `annotation_df`.", call. = FALSE)
    }
    out$strand <- as.character(annotation_df[[strand_col]])
  } else {
    out$strand <- "."
  }

  if (!is.null(attributes_col)) {
    if (!(attributes_col %in% names(annotation_df))) {
      stop("`attributes_col` is not present in `annotation_df`.", call. = FALSE)
    }
    out$attributes <- as.character(annotation_df[[attributes_col]])
  }

  if (any(is.na(out$seqname) | out$seqname == "")) {
    stop("`seqname` must not contain missing or empty values.", call. = FALSE)
  }
  if (any(is.na(out$feature_type) | out$feature_type == "")) {
    stop("`feature_type` must not contain missing or empty values.", call. = FALSE)
  }
  if (any(is.na(out$feature_id) | out$feature_id == "")) {
    stop("`feature_id` must not contain missing or empty values.", call. = FALSE)
  }
  if (any(is.na(out$start) | is.na(out$end) | out$start <= 0 | out$end < out$start)) {
    stop("`start`/`end` must be positive integers with `end >= start`.", call. = FALSE)
  }

  valid_strand <- c("+", "-", ".")
  out$strand[is.na(out$strand) | out$strand == ""] <- "."
  if (any(!(out$strand %in% valid_strand))) {
    stop("`strand` must contain only '+', '-', or '.'.", call. = FALSE)
  }

  out
}

#' Read GFF/GTF-style annotation table into canonical feature schema
#'
#' Parses standard 9-column GFF/GTF files and returns CoMMA's canonical feature
#' table. `feature_id` is derived from `ID=` when present, then from `Name=`, and
#' finally falls back to a deterministic `feature_type:start-end` identifier.
#'
#' @param path Path to a GFF/GTF file.
#' @param comment.char Character indicating comment lines (default `"#"`).
#'
#' @return A canonical feature table with required columns and optional
#'   `attributes`.
#' @export
read_annotation_gff <- function(path, comment.char = "#") {
  gff <- utils::read.delim(
    file = path,
    header = FALSE,
    sep = "\t",
    comment.char = comment.char,
    stringsAsFactors = FALSE,
    quote = ""
  )

  if (ncol(gff) < 9) {
    stop("GFF input must contain at least 9 tab-delimited columns.", call. = FALSE)
  }

  gff <- gff[, 1:9]
  names(gff) <- c(
    "seqname", "source", "feature_type", "start", "end",
    "score", "strand", "phase", "attributes"
  )

  extract_attr <- function(attr, key) {
    pattern <- paste0("(?:^|;)", key, "=([^;]+)")
    matches <- regexec(pattern, attr)
    out <- regmatches(attr, matches)
    value <- rep(NA_character_, length(attr))
    has_match <- lengths(out) > 0
    value[has_match] <- vapply(out[has_match], function(x) x[2], character(1))
    value
  }

  feature_id <- extract_attr(gff$attributes, "ID")
  name_id <- extract_attr(gff$attributes, "Name")
  fill_idx <- is.na(feature_id) | feature_id == ""
  feature_id[fill_idx] <- name_id[fill_idx]
  fill_idx <- is.na(feature_id) | feature_id == ""
  feature_id[fill_idx] <- paste0(
    gff$feature_type[fill_idx], ":", gff$start[fill_idx], "-", gff$end[fill_idx]
  )

  normalize_annotation_table(
    annotation_df = data.frame(gff, feature_id = feature_id, stringsAsFactors = FALSE),
    feature_type_col = "feature_type",
    feature_id_col = "feature_id",
    seqname_col = "seqname",
    start_col = "start",
    end_col = "end",
    strand_col = "strand",
    attributes_col = "attributes"
  )
}

#' Read BED-style annotation table into canonical feature schema
#'
#' Parses common BED tables. BED start is converted from 0-based to 1-based
#' coordinates. If no name column exists, deterministic feature IDs are
#' generated as `feature_type:start-end`.
#'
#' @param path Path to a BED file.
#' @param feature_type Default feature type to assign when a BED type column is
#'   unavailable.
#' @param comment.char Character indicating comment lines (default `"#"`).
#'
#' @return A canonical feature table with required columns and optional
#'   `attributes`.
#' @export
read_annotation_bed <- function(path,
                                feature_type = "feature",
                                comment.char = "#") {
  bed <- utils::read.delim(
    file = path,
    header = FALSE,
    sep = "\t",
    comment.char = comment.char,
    stringsAsFactors = FALSE,
    quote = ""
  )

  if (ncol(bed) < 3) {
    stop("BED input must contain at least 3 tab-delimited columns.", call. = FALSE)
  }

  seqname <- as.character(bed[[1]])
  start <- as.integer(bed[[2]]) + 1L
  end <- as.integer(bed[[3]])

  feature_id <- if (ncol(bed) >= 4) as.character(bed[[4]]) else NA_character_
  strand <- if (ncol(bed) >= 6) as.character(bed[[6]]) else "."

  feature_type_vec <- rep(feature_type, nrow(bed))

  missing_id <- is.na(feature_id) | feature_id == ""
  feature_id[missing_id] <- paste0(feature_type_vec[missing_id], ":", start[missing_id], "-", end[missing_id])

  attributes <- NULL
  if (ncol(bed) > 6) {
    extra <- bed[, 7:ncol(bed), drop = FALSE]
    attributes <- apply(extra, 1, function(x) paste(x, collapse = ";"))
  }

  parsed <- data.frame(
    feature_type = feature_type_vec,
    feature_id = feature_id,
    seqname = seqname,
    start = start,
    end = end,
    strand = strand,
    stringsAsFactors = FALSE
  )

  if (!is.null(attributes)) {
    parsed$attributes <- attributes
    normalize_annotation_table(
      parsed,
      feature_type_col = "feature_type",
      feature_id_col = "feature_id",
      seqname_col = "seqname",
      start_col = "start",
      end_col = "end",
      strand_col = "strand",
      attributes_col = "attributes"
    )
  } else {
    normalize_annotation_table(
      parsed,
      feature_type_col = "feature_type",
      feature_id_col = "feature_id",
      seqname_col = "seqname",
      start_col = "start",
      end_col = "end",
      strand_col = "strand"
    )
  }
}

#' Annotate methylation sites by overlap with canonical features
#'
#' Performs base-resolution overlap between methylation site calls and canonical
#' features. Many-to-one relationships are preserved: each site can appear in
#' multiple output rows if it overlaps multiple features.
#'
#' @param site_table Site table with at least `seqname` and position columns.
#' @param feature_table Canonical feature table from
#'   [normalize_annotation_table()], [read_annotation_gff()], or
#'   [read_annotation_bed()].
#' @param site_pos_col Column name in `site_table` containing 1-based positions.
#' @param site_seqname_col Column name in `site_table` containing sequence names.
#' @param site_strand_col Optional column name for site strand. Required when
#'   `match_strand = TRUE`.
#' @param match_strand If `TRUE`, only overlaps with matching strand are kept.
#'   Features with `strand == "."` are treated as strand-agnostic.
#' @param include_unannotated If `TRUE`, emit one row for unmatched sites with
#'   `NA` feature fields.
#'
#' @return A long `data.frame` with original site columns and appended feature
#'   columns.
#' @export
annotate_sites_with_features <- function(site_table,
                                         feature_table,
                                         site_pos_col = "pos",
                                         site_seqname_col = "seqname",
                                         site_strand_col = NULL,
                                         match_strand = FALSE,
                                         include_unannotated = TRUE) {
  if (!is.data.frame(site_table)) {
    stop("`site_table` must be a data.frame.", call. = FALSE)
  }
  if (!is.data.frame(feature_table)) {
    stop("`feature_table` must be a data.frame.", call. = FALSE)
  }

  required_feature <- c("feature_type", "feature_id", "seqname", "start", "end", "strand")
  missing_feature <- setdiff(required_feature, names(feature_table))
  if (length(missing_feature) > 0) {
    stop(
      paste0("`feature_table` is missing required columns: ", paste(missing_feature, collapse = ", "), "."),
      call. = FALSE
    )
  }

  if (!(site_pos_col %in% names(site_table))) {
    stop("`site_pos_col` is not present in `site_table`.", call. = FALSE)
  }
  if (!(site_seqname_col %in% names(site_table))) {
    stop("`site_seqname_col` is not present in `site_table`.", call. = FALSE)
  }
  if (match_strand && (is.null(site_strand_col) || !(site_strand_col %in% names(site_table)))) {
    stop("`site_strand_col` must be provided and present when `match_strand = TRUE`.", call. = FALSE)
  }

  site_idx <- seq_len(nrow(site_table))
  out_rows <- vector("list", length = 0)

  for (i in site_idx) {
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
