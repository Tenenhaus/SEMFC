# correction of the RGCCA estimates


# source("R/utils/cov2.R")


# correction = function(A, a, bias = FALSE){
#   S = cov2(A, bias = bias)
#   diag(S) = 0
#   d = drop(sqrt((t(a)%*%S%*%a)/(1-sum(a^4))))
#   return(list(d = d, n = 0))
# }


#' Correct Block Weight Estimates in Component-Based SEM
#'
#' Computes a correction factor for block weights obtained from RGCCA-like algorithms.
#' The correction differs depending on whether the measurement model is reflective
#' or formative.
#'
#' @param A Numeric matrix of indicators for the block (n × p).
#' @param a Numeric vector of block weights (length p).
#' @param mode Character string indicating the measurement mode: "reflective" or "formative".
#'   Default is "reflective".
#' @param bias Logical indicating whether to use biased covariance estimator (division by n)
#'   or unbiased (division by n-1). Default is `FALSE` (unbiased).
#'
#' @return List containing:
#'   \item{d}{Numeric correction factor to apply to the block weights.}
#'   \item{n}{Integer count of negative covariances (reflective) or zero eigenvalues (formative).}
#'
#' @details
#' ## Reflective mode
#' For reflective blocks, the correction factor is computed as:
#' \deqn{d = \sqrt{\frac{\sum_{i,j: (W \odot S)_{ij} > 0} (W \odot S)_{ij}}{\sum_{i,j: (W \odot S)_{ij} > 0} W_{ij}^2}}}
#' where \eqn{W = aa'} (with diagonal set to 0) and \eqn{S} is the sample covariance matrix.
#' If all elements are negative, \eqn{d = 1}.
#'
#' ## Formative mode
#' For formative blocks, the correction factor is:
#' \deqn{d = \frac{1}{\sqrt{a'S^{-1}a}}}
#' where \eqn{S} is the sample covariance matrix of indicators.
#'
#' @keywords internal
correction = function(A, a, mode = "reflective", bias = FALSE){
  
  if(mode == "reflective"){
    S = cov2(A, bias = bias)
    W = a%*%t(a);diag(W) = 0
    WS = W*S
    d = ifelse(any(WS>0), sqrt(sum(WS[WS>0])/sum(W[WS>0]^2)), 1)
    #d =  sqrt(sum(WS)/sum(W^2))
    n = length(which(WS<0))
  }
  
  if(mode == "formative"){
    #S = cov2(A, bias = bias)
    #eig = eigen(S)
    #Sinv_sqrt  = eig$vectors %*% diag(eig$values^(-1/2)) %*% t(eig$vectors)
    #d = 1/norm(Sinv_sqrt%*%a, type = "2")
    #n = sum(eig$values == 0)
    d = 1/drop(sqrt(t(a)%*%solve(cov2(A, bias = bias))%*%a))
    n = 0
  }
  
  return(list(d = d, n = n))
}
  
  
