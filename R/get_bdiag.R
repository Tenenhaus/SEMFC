#' Construct Block Diagonal Covariance Structures for Formative and Reflective Blocks
#'
#' This function creates a list of block diagonal covariance matrix for both formative and reflective blocks
#' using the parameter vector `x`, block sizes, and the starting index of covariance parameters.
#'
#' @param x A numeric vector of parameters from which covariance or variance values are extracted.
#' @param mode A character vector indicating the mode ("formative" or "reflective") for each block.
#' @param block_sizes A numeric vector where each element specifies the number of indicators
#'   of each block.
#' @param initial_start_index_cov An integer indicating the starting index in `x` for covariance/variance
#'   parameters.
#'
#' @return A list of block matrices, where each block corresponds to either:
#'   - A covariance matrix for formative blocks.
#'   - A diagonal matrix for reflective blocks.
#'
#' @examples
#' \dontrun{
#' set.seed(27)
#' x <- rnorm(61)
#' mode <- c("formative", "formative", "formative", "formative", "reflective", "reflective")
#' block_sizes <- rep(3, 6)
#' initial_start_index_cov <- 32
#' get_bdiag(x, mode, block_sizes, initial_start_index_cov)}
#' @keywords internal

get_bdiag <- function(x, mode, block_sizes, initial_start_index_cov) {


  # number of blocks
  J <- length(block_sizes)
  # list of lengths of the upper values in the cov matrix for the composite block i or
  # of the diagonal values for formative for each block
  lengths_values_cov <- block_sizes
  lengths_values_cov[mode == "formative"] <- (block_sizes[mode == "formative"]^2 + block_sizes[mode == "formative"]) / 2
  # number of parameters for covariance
  total_cov_parameter <- sum(lengths_values_cov)
  end_endex_cov <- initial_start_index_cov + total_cov_parameter - 1

  # part of the vector corresponding to covariance blocks
  extracted_parameters_cov <- x[initial_start_index_cov:end_endex_cov]
  # list of parameters corresponding to each covariance bloc
  list_cov <- split(extracted_parameters_cov,
                    rep(seq_along(lengths_values_cov), lengths_values_cov))
  names(list_cov) <- names(lengths_values_cov)

  # Building list of formative matrices
  S_composites <- lapply(list_cov[mode == "formative"], build_formative_S_diag)

  # Building list of reflectives matrices
  reflective_blocks <- lapply(list_cov[mode == "reflective"], diag)

  BDIAG <- vector("list", J)
  BDIAG[mode == "formative"] <- S_composites
  BDIAG[mode == "reflective"] <- reflective_blocks
  names(BDIAG) <- names(lengths_values_cov)

  return(BDIAG)

}