


#' Construct Correlation Matrix from Parameter Vector
#'
#' This function constructs a symmetric correlation matrix for a set of latent variables
#' based on the parameter vector `x` and the starting index for correlation-related parameters.
#'
#' @param x A numeric vector of parameters. Correlation-related values are extracted starting
#'   from the index specified by `start_index`.
#' @param latent_variables A vector of names or identifiers for the latent variables, used
#'   to determine the dimensions of the correlation matrix.
#' @param start_index An integer specifying the index in `x` where the correlation parameters
#'   begin.
#'
#' @return A symmetric correlation matrix (dim x dim), where `dim` is the number of latent variables.
#'
#' @examples
#' latent_variables <- c("LV1", "LV2", "LV3")
#' x <- c(0.5, 0.3, 0.2, 0.5, 0.4) # upper triangle values for a 3x3 matrix
#' start_index <- 1
#' get_correlation_coeff(x, latent_variables, start_index)


get_correlation_coeff <- function(x, latent_variables, start_index) {


  # dim  matrix (number of exogeneous or endogeneous latent variables)
  dim <- length(latent_variables)
  #length upper values vector from this symetric matrix:
  length_correl <- dim * (dim - 1) / 2
  end_index <- start_index + length_correl - 1
  upper_values <- x[start_index:end_index]

  # Building of the correlation matrix
  P_matrix <- diag(1, dim, dim)
  if (dim > 1) {
    P_matrix[upper.tri(P_matrix)] <- upper_values
    P_matrix[lower.tri(P_matrix)] <- t(P_matrix)[lower.tri(P_matrix)]
  }

  return(P_matrix)
}
