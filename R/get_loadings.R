#' Extract Loadings from Parameter Vector
#'
#' This function retrieves the loadings for each block of observed variables
#' based on the parameter vector `x` and the input block sizes.
#'
#' @param x Numeric vector containing all model parameters (loadings, correlations,
#'   path coefficients, and variance/covariance parameters).
#' @param block_sizes A numeric vector where each element represents the number
#'   of observed variables in the corresponding block of latent variables.
#'
#' @return A list of numeric vectors, where each vector contains the loadings
#'   for the corresponding block of manifest variables.
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' block_sizes <- c(3, 3, 3, 3, 3, 3)
#' x  <- rnorm(61)
#' get_loadings(x, block_sizes)}
#' @keywords internal
get_loadings <- function(x, block_sizes) {


  # total number of emergent variables
  y <- sum(block_sizes)
  loadings <- setNames(
    lapply(seq_along(block_sizes), function(b) {
      setNames(
        split(x[1:y], rep(seq_along(block_sizes), block_sizes))[[b]],
        paste0(names(block_sizes[b]), 1:block_sizes[b])
      )
    }),
    names(block_sizes)
  )


  return(loadings)

}