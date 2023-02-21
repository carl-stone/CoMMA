#' Symmetrical Logarithm
#'
#' @param x a numeric or complex vector
#' @param base a positive or complex number, passed to base [log()]. Defaults
#'  to *e*`=exp(1)`.
#' @param C a numeric.
#'
#' @return A vector of the same length as `x` containing the transformed values.
#'  Unlike `log`, all values `x` are defined for `symlog(x)`.
#' @export
#'
#' @details
#' `y = sign(x) * log(1 + (abs(x) / 10^C))`
#' where `C` is a constant that affects how the slope of the function through
#' the origin. With `C = 0` the slope at the origin is 1.
#'
#' While the output values between `log` and `symlog` are not identical,
#' `symlog` has unique behavior that can be useful in specific situations.
#'
#' @references
#' [J Beau W Webber 2013 Meas. Sci. Technol. 24 027001](https://doi.org/10.1088/0957-0233/24/2/027001)
#'
#' @examples
#' log(-10:10)
#' # Returns NaN or -Inf for -10:0
#'
#' symlog(-10:10)
#' # All values can be computed.
#'
#' # The difference between symlog() and log() is max at x = 1 (delta = 0.693)
#' # and quickly approaches 0.
#' plot(1:100, symlog(1:100) - log(1:100))
symlog <- function(x, base = exp(1), C = 0) {
  sign(x) * log(1 + (abs(x) / 10^C), base = base)
}

#' Convert Methylation B-values to M-values
#'
#' @param B a numeric vector of Beta-values.
#' @param alpha a double. Used to prevent division by zero and should be an
#'  order of magnitude smaller than your smallest beta value.
#' @param sym a boolean. If `TRUE`, use the symmetrical logarithm [symlog()]. If
#' `FALSE`, use base function [log()].
#' @param base a positive numeric. The base with respect to which logarithms
#'  are computed. Ultimately passed to [log()]. Don't change this.
#'
#'
#' @return a numeric vector of M-values
#' @export
#'
#' @examples
#' # Sample data to make a mock methylation sample.
#' data <- sample(all_samples_long$beta, size = 10000, replace = TRUE)
#'
#' # Use default log2 when all beta values fall between 0 and 1.
#' beta <- methylBtoM(data)
#'
#' # Use symmetrical log2 when values fall outside 0 and 1.
#' beta <- methylBtoM(data, sym = TRUE)
methylBtoM <- function(B, alpha = 0.001, sym = FALSE, base = 2) {
  if (sym == FALSE) {
    log((B + alpha) / (1 - B + alpha), base = base)
  }
  else if (sym == TRUE) {
    symlog((B) / (1 - B), base = base)
  }
}

#' Title Convert Methylation M-values to B-values
#'
#' @param M a numeric vector
#'
#' @return a numeric vector
#' @export
#'
#' @examples
#' data <- sample(all_samples_long$beta, size = 10000, replace = TRUE)
#' mvalues <- methylBtoM(data)
#' bvalues <- methylMtoB(mvalues)
methylMtoB <- function(M) {
  2^M / (1 + 2^M)
}
