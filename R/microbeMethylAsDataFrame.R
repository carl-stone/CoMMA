microbeMethylAsDataFrame <- function(mm) {
  mm <- as.data.frame(mm)
  names(mm)[names(mm) %in% c("group", "group_name")] <- c("sample", "sample_name")
}
