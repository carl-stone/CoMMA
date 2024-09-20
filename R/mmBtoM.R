mmBtoM <- function(mm) {
  min_B <- assays(mm)$B[assays(mm)$B > 0] |> min(na.rm = T)
  assays(mm)$M <- methylBtoM(assays(mm)$B, alpha = min_B)

  return(mm)
}
