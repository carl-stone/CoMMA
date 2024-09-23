mmFilter <- function(mm, expr, na.action = c("omit", "include", "fail")) {
  # Match na.action argument
  na.action <- match.arg(na.action)

  expr <- substitute(expr)

  rr_dat <- as.data.frame(unlist(rowRanges(mm)), row.names = NULL)

  logical_vec <- eval(expr, envir = rr_dat)

  # Handle NAs in logical_vec according to na.action parameter
  if (any(is.na(logical_vec))) {
    if (na.action == "omit") {
      logical_vec[is.na(logical_vec)] <- FALSE
    } else if (na.action == "include") {
      logical_vec[is.na(logical_vec)] <- TRUE
    } else if (na.action == "fail") {
      stop("logical subscript contains NAs")
    }
  }

  q <- rr_dat[logical_vec, ] |>
    as_granges()

  mm_filtered <- subsetByOverlaps(mm, q)

  return(mm_filtered)

}

