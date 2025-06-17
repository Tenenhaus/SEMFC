#' Extract Loadings from Parameter Vector
#'
#' This function retrieves the loadings for each block of observed variables
#' (manifest variables) based on the parameter vector `x` and the input block sizes.
#'
#' @param x A numeric vector of parameters. The first `y` values in `x` represent
#'   the loadings to be extracted, where `y` is the total number of manifest
#'   variables across all blocks.
#' @param block_sizes A numeric vector where each element represents the number
#'   of observed variables (manifest variables) in the corresponding block of
#'   latent variables.
#'
#' @return A list of numeric vectors, where each vector contains the loadings
#'   for the corresponding block of manifest variables.
#'
#' @examples
#' block_sizes <- c(3, 2)
#' x <- c(0.2, 0.4, 0.6, 0.1, 0.3)
#' get_loadings(x, block_sizes)

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