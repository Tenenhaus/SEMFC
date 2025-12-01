
#' Standardized Root Mean Square Residual (SRMR) for Structural Equation Models
#'
#' This function calculates the SRMR, a goodness-of-fit measure for structural equation models,
#' by comparing the empirical covariance matrix with the model-implied covariance matrix.
#'
#' @param S A covariance matrix representing the empirical data.
#' @param Sigma A covariance matrix representing the model-implied data.
#'
#' @return A numeric value representing the SRMR, which is the square root of the mean
#'   of the squared residuals between the lower triangular elements (including the diagonal)
#'   of the standardized empirical and model-implied correlation matrices.
#'
#' @details The function first standardizes the covariance matrices to correlation matrices
#'   using \code{cov2cor}. It then computes the residuals as the difference between the
#'   empirical and model-implied correlation matrices. The SRMR is calculated as the square
#'   root of the mean of the squared residuals from the lower triangular part of the matrix.
#'
#' @examples
#' \dontrun{
#' S <- matrix(c(1, 0.5, 0.5, 1), ncol = 2)
#' Sigma <- matrix(c(1, 0.4, 0.4, 1), ncol = 2)
#' srmr <- srmrsem(S, Sigma)
#' print(srmr)
#' }
#'
#' @export
srmrsem <- function(S, Sigma) {
  R_emp <- cov2cor(S)
  R_mod <- cov2cor(Sigma)
  resids <- R_emp - R_mod
  SRMR <- sqrt(mean(resids[lower.tri(resids, diag=TRUE)]^2))
  return(SRMR)
}

