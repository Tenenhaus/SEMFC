


source("R/ml_sem/lvm_ml/build_formative_S_diag.R")


#' Construct Block Diagonal Covariance Structures for Formative and Reflective Blocks
#'
#' This function creates a block diagonal covariance matrix for both formative and reflective blocks
#' using the parameter vector `x`, block sizes, and the starting index of covariance parameters.
#'
#' @param x A numeric vector of parameters from which covariance or variance values are extracted.
#' @param mode A character vector indicating the mode ("formative" or "reflective") for each block.
#' @param block_sizes A numeric vector where each element specifies the size (number of variables)
#'   of each block.
#' @param loadings A list of numeric vectors representing loadings for each block. Each vector corresponds
#'   to the respective block size.
#' @param initial_start_index_cov An integer indicating the starting index in `x` for covariance/variance
#'   parameters.
#'
#' @return A list of block diagonal matrices, where each block corresponds to either:
#'   - A covariance matrix for formative blocks.
#'   - A diagonal matrix for reflective blocks.
#'
#' @examples
#' x <- c(1, 0.3, 1, 0.2, 0.4, 1, 2, 0.5, 3)
#' mode <- c("formative", "reflective")
#' block_sizes <- c(2, 3)
#' loadings <- list(c(0.8, 0.6), c(0.7, 0.5, 0.4))
#' initial_start_index_cov <- 1
#' get_bdiag_bis(x, mode, block_sizes, loadings, initial_start_index_cov)


get_bdiag_bis <- function(x, mode, block_sizes, initial_start_index_cov) {


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

  # Building list of formative matrices
  S_composites <- lapply(list_cov[mode == "formative"], build_formative_S_diag)

  # Building list of reflectives matrices
  reflective_blocks <- lapply(list_cov[mode == "reflective"], diag)

  BDIAG <- vector("list", J)
  BDIAG[mode == "formative"] <- S_composites
  BDIAG[mode == "reflective"] <- reflective_blocks

  return(BDIAG)

}