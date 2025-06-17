
#' Construct Symmetric Formative S Matrix
#'
#' This function builds a symmetric covariance matrix (S_composite_i) for a formative block,
#' using the upper triangular values provided in `upper_values_composite_i`.
#'
#' @param upper_values_composite_i A numeric vector containing the values for the upper
#'   triangle (including the diagonal) of the symmetric matrix.
#'
#' @return A symmetric formative covariance matrix of dimensions `dim_block x dim_block`.
#'
#' @examples
#' upper_values <- c(1, 0.3, 0.5, 1, 0.2, 1)
#' dim_block <- 3
#' build_formative_S_diag(upper_values, dim_block)



build_formative_S_diag <- function(upper_values_composite_i, dim_block) {


  # dim of the matrix extracted from the length of the upper values
  dim_block <- (-1 + sqrt(1 + 8 * length(upper_values_composite_i))) / 2
  S_composite_i <- matrix(0, dim_block, dim_block)
  S_composite_i[upper.tri(S_composite_i, diag = TRUE)] <- upper_values_composite_i
  S_composite_i[lower.tri(S_composite_i)] <- t(S_composite_i)[lower.tri(S_composite_i)]

  return(S_composite_i)
}
