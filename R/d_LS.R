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
#' S <- diag(2)
#' SIGMA <- diag(2) * 0.5
#' discrepancy <- d_LS(S, SIGMA)
#'
#' @export
d_LS <- function(A, B){
  A <- as.matrix(A)
  B <- as.matrix(B)
  d_LS <- sqrt(0.5*sum(diag((B-A)%*%(B-A))))
  return(d_LS = d_LS)
}