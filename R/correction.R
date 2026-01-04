#' Correct Loading Estimates in Component-Based SEM
#'
#' Computes the loading norm estimator giving the normalised loadings.
#' The correction differs depending on whether the measurement model is reflective
#' or formative.
#'
#' @param A Numeric matrix of indicators for the block (n × p).
#' @param a Numeric vector of normalised loadings (length p).
#' @param mode Character string indicating the measurement mode: "reflective" or "formative".
#'   Default is "reflective".
#' @param bias Logical indicating whether to use biased covariance estimator (division by n)
#'   or unbiased (division by n-1). Default is `FALSE` (unbiased).
#'
#' @return List containing:
#'   \item{d}{Numeric correction factor to apply to the loadings.}
#'   \item{n}{Integer count of negative covariances (reflective) or zero eigenvalues (formative).}
#'
#' @details
#' ## Reflective mode
#' For reflective blocks, the correction factor is computed as:
#' \deqn{d = \sqrt{
#'   \frac{\mathbf{a}_j^t (\mathbf{S}_{jj} - \text{diag}(\mathbf{S}_{jj})) \mathbf{a}_j}
#'   {1 - \mathbf{a}_j^t \text{diag}(\mathbf{a}_j \mathbf{a}_j^t) \mathbf{a}_j}
#' }}
#' where \eqn{S} is the sample covariance matrix.
#' If all elements are negative, \eqn{d = 1}.
#'
#' ## Formative mode
#' For formative blocks, the correction factor is:
#' \deqn{d = \frac{1}{\sqrt{a'S^{-1}a}}}
#' where \eqn{S} is the sample covariance matrix of indicators.
#'
#' @keywords internal
correction <- function(A, a, mode = "reflective", bias = FALSE){
  
  if(mode == "reflective"){
    S <- cov2(A, bias = bias)
    W <- a%*%t(a);diag(W) <- 0
    WS <- W*S
    d <- ifelse(any(WS>0), sqrt(sum(WS[WS>0])/sum(W[WS>0]^2)), 1)
    n <- length(which(WS<0))
  }
  
  if(mode == "formative"){
    #S = cov2(A, bias = bias)
    #eig = eigen(S)
    #Sinv_sqrt  = eig$vectors %*% diag(eig$values^(-1/2)) %*% t(eig$vectors)
    #d = 1/norm(Sinv_sqrt%*%a, type = "2")
    #n = sum(eig$values == 0)
    d <- 1/drop(sqrt(t(a)%*%solve(cov2(A, bias = bias))%*%a))
    n <- 0
  }
  
  return(list(d = d, n = n))
}
  
  
