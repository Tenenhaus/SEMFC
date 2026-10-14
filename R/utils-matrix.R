
duplication_matrix <- function(P) {

  n_vech <- P * (P + 1) / 2
  mat_index <- matrix(0, nrow = P, ncol = P)
  mat_index[lower.tri(mat_index, diag = TRUE)] <- 1:n_vech
  mat_index[upper.tri(mat_index)] <- t(mat_index)[upper.tri(mat_index)]
  index_col <- as.vector(mat_index)
  D_P <- sparseMatrix(
    i = 1:(P^2),
    j = index_col,
    x = 1,
    dims = c(P^2, n_vech)
  )

  return(D_P)
}


correlation_duplication_matrix <- function(n, D_n = duplication_matrix(n)) {

  # we remove the columns corresponding to the diagonal elements in vech() to get the correlation duplication matrix
  # Direct computation of the indices of the diagonal elements in vech() for a matrix of size n x n
  j <- 1:n
  indices_diag <- 1 + (j - 1) * n - (j - 1) * (j - 2) / 2
  return(D_n[, -indices_diag, drop = FALSE])
}



commutation_matrix <- function(m, n = m) {
  mn <- m * n
  index_col <- as.vector(t(matrix(1:mn, nrow = m, ncol = n)))
  K_mn <- sparseMatrix(
    i = 1:mn,
    j = index_col,
    x = 1,
    dims = c(mn, mn)
  )
  return(K_mn)
}

diagonal_extraction_matrix <- function(m) {
  L_D <- sparseMatrix(
    i = 1:m,
    j = seq(from = 1, to = m^2, by = m + 1),
    x = 1,
    dims = c(m, m^2)
  )
  return(L_D)
}

elimination_matrix<- function(m) {

  nb_vech <- m * (m + 1) / 2
  mat_index <- matrix(1:(m^2), nrow = m, ncol = m)
  idx_keep <- mat_index[lower.tri(mat_index, diag = TRUE)]
  L_m <- sparseMatrix(
    i = 1:nb_vech,
    j = idx_keep,
    x = 1,
    dims = c(nb_vech, m^2)
  )
  return(L_m)
}



#' Least Squares Discrepancy Function
#'
#' Computes the Euclidean distance between two covariance matrices using the
#' least squares discrepancy function.
#'
#' @param A A matrix, typically the observed covariance matrix
#' @param B A matrix, typically the model-implied covariance matrix
#'
#' @return The least squares discrepancy value representing the
#'   Euclidean distance between matrices A and B
#'
#' @details The LS discrepancy function is defined as:
#'   \deqn{d_{LS} = \sqrt{frac{1}{2} \times \text{tr}[(B-A)^2]}}
#'   where tr denotes the trace operator. This measures the overall fit
#'   between the observed and model-implied covariance matrices.
#'
#' @examples
#' \dontrun{
#' S <- diag(2)
#' SIGMA <- diag(2) * 0.5
#' discrepancy <- d_LS(S, SIGMA)
#'}
#'
#' @keywords internal
d_LS <- function(A, B){
  A <- as.matrix(A)
  B <- as.matrix(B)
  d_LS <- sqrt(0.5*sum(diag((B-A)%*%(B-A))))
  return(d_LS = d_LS)
}



#' Compute Reliability for Formative Blocks
#'
#' This function calculates the reliability of formative measurement blocks
#' using the specified metric. Currently, only the "Dillon" metric is supported.
#'
#' @param metric A character string specifying the reliability metric to use.
#'   Currently, only "Dillon" is implemented.
#' @param lambdas A list of factor loadings for each block.
#' @param residual_variances A list of residual variances for each block.
#'
#' @return A named numeric vector containing the reliability values for each block.
#'
#' @details The Dillon reliability metric is calculated as the squared sum of
#'   loadings divided by the sum of the squared loadings and the residual variances.
#'
#' @examples
#' \dontrun{
#' lambdas <- list(block1 = c(0.8, 0.7), block2 = c(0.9, 0.6, 0.5))
#' residuals <- list(block1 = c(0.2, 0.3), block2 = c(0.1, 0.4, 0.5))
#' reliability_values <- reliability(
#'  metric = "Dillon",
#'  lambdas = lambdas,
#'  residual_variances = residuals
#' )
#' print(reliability_values)
#' }
#' @keywords internal

# pour bloc formatif
reliability <- function (metric, lambdas, residual_variances){
  if (metric=='Dillon'){
    # compute sum of loadings and residual variances for each block
    sum_loadings <- sapply(lambdas, function(x) sum(as.numeric(x)))
    sum_residuals <- sapply(residual_variances, function(x) sum(as.numeric(x)))

    #remove formative blocks (blocks without residual variances)
    sum_loadings <- sum_loadings[names(sum_loadings) %in% names(sum_residuals)]

    # compute numerator and denominator for Dillon's rho
    numerator <- sum_loadings^2
    denominator <- numerator + sum_residuals

    # compute reliability results
    results <- numerator / denominator

    return(results)
  }
}

