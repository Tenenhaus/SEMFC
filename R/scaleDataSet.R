

#' Scale Data Set Based on Covariance Matrices
#'
#' This function transforms a data set using the method proposed by Yuan & Hayashi (2003),
#' which involves scaling the data based on the covariance matrix of the data and a target covariance matrix.
#'
#' @param data A data frame or matrix representing the data set to be scaled.
#' @param Sigma A covariance matrix to which the data will be scaled.
#' @param bias A logical value indicating whether to use a biased covariance matrix estimate.
#'   Defaults to \code{FALSE}.
#'
#' @return A matrix representing the scaled data set, with the same column names as the input data.
#'
#' @details The function computes the covariance matrix of the data, performs singular value decomposition (SVD)
#'   on both the data covariance matrix and the target covariance matrix, and applies the transformation:
#'   \deqn{ScaledData = Data \times S^{-1/2} \times \Sigma^{1/2}}
#'   where \eqn{S^{-1/2}} is the inverse square root of the data covariance matrix and \eqn{\Sigma^{1/2}}
#'   is the square root of the target covariance matrix.
#'
#' @examples
#' \dontrun{
#' data <- matrix(rnorm(100), ncol = 5)
#' Sigma <- diag(5)
#' scaled_data <- scaleDataSet(data, Sigma)
#' }
#'
#' @export

scaleDataSet <- function(data, Sigma, bias = FALSE){
  S <- cov2(data, bias = bias)
  # singular value decomposition
  S1 = svd(S); S2 = svd(Sigma)
  d1 = S1$d; d2 = S2$d
  u1 = S1$u; u2 = S2$u
  v1 = S1$v; v2 = S2$v
  S_half <- u1%*%(diag(d1^(-1/2)))%*%t(v1) #S^(-1/2)
  Sigma_half <- u2%*%(diag(d2^(1/2)))%*%t(v2) # Sigma^(1/2)
  ScaleData <- as.matrix(data)%*%S_half%*%Sigma_half
  colnames(ScaleData) <- colnames(data)
  return(ScaleData = ScaleData)
}
