#' cov2() is similar to cov() but has an additional argument. The denominator \eqn{n} (bias = TRUE)
#' can be used (instead of \eqn{n-1}) to give a biased estimator of the (co)variance.
#' @param x A numeric vector, matrix or data.frame.
#' @param y A numeric vector, matrix or data.frame.
#' @param bias A logical value. If bias = TRUE, \eqn{n} is used to give a biased estimator of the (co)variance.
#' If bias = FALSE, \eqn{n-1} is used (default: TRUE).
#' @return \item{C}{Estimation of the variance (resp. covariance) of x (resp. x and y).}
#' @title Variance and Covariance (Matrices)
#' @importFrom stats cov
#' @keywords internal
cov2 <- function (x, y = NULL, bias = TRUE)
{
  n <- NROW(x)

  if (is.null(y)) {
    x <- as.matrix(x)
    if(bias){
      C <- ((n - 1)/n) * cov(x, use = "pairwise.complete.obs")
    }
    else{
      C <- cov(x, use = "pairwise.complete.obs")
    }
  }

  else{
    if(bias){
      C <- ((n - 1)/n) * cov(x, y, use = "pairwise.complete.obs")
    }
    else{
      C <- cov(x, y, use = "pairwise.complete.obs")
    }
  }
return(C)
}




#' Standardization (to zero means and unit variances) of matrix-like objects.
#' @param A A numeric matrix.
#' @param center A logical value. If center = TRUE, each column is translated to have zero mean (default: TRUE).
#' @param scale A logical value. If scale = TRUE, each column is transformed to have unit variance (default = TRUE).
#' @param bias Logical value for biaised (\eqn{1/n}) or unbiaised (\eqn{1/(n-1)}) estimator of the var/cov (default = TRUE).
#' @return \item{A}{The centered and/or scaled matrix. The centering and scaling values (if any) are returned as attributes "scaled:center" and "scaled:scale".}
#' @title Scaling and Centering of Matrix-like Objects
#' @keywords internal

scale2 <- function (A, center = TRUE, scale = TRUE, bias = TRUE)
{
  if (center == TRUE & scale == TRUE) {
    A <- scale(A, center = TRUE, scale = FALSE)
    std <- sqrt(apply(A, 2, function(x) cov2(x, bias = bias)))
    if (any(std==0)) {
      sprintf("there were %d constant variables",sum(std==0))
      std[std==0] <- 1
    }
    A <- A/matrix(rep(std, NROW(A)), NROW(A), NCOL(A), byrow = TRUE)
    attr(A, "scaled:scale") <- std
    return(A)
  }
  if (center == TRUE & scale == FALSE) {
    A <- scale(A, center = TRUE, scale = FALSE)
    return(A)
  }
  if (center == FALSE & scale == TRUE) {
    std <- apply(A, 2, function(x) cov2(x, bias = bias))
    A <- A/matrix(rep(std, NROW(A)), NROW(A), NCOL(A), byrow = TRUE)
    attr(A, "scaled:scale") <- std
    return(A)
  }
}


pm <- function(M1, M2, na.rm = TRUE) {
  if (na.rm) {
    M1[is.na(M1)] <- 0
    M2[is.na(M2)] <- 0
  }
  return(M1 %*% M2)
}




tau_estimate <- function(x, na.rm = TRUE) {
  p <- NCOL(x)
  if (p == 1) {
    return(1)
  }
  n <- NROW(x)
  corm <- cor(x, use = "pairwise.complete.obs")
  xs <- scale(x, center = TRUE, scale = TRUE)
  v <- (n / ((n - 1) ^ 3)) * (pm(t(xs ^ 2), xs ^ 2, na.rm = na.rm) - 1 / n * pm(t(xs), xs, na.rm = na.rm) ^ 2)
  diag(v) <- 0
  I <- diag(NCOL(x))
  d <- (corm - I) ^ 2
  tau <- sum(v) / sum(d)
  tau <- max(min(tau, 1), 0)
  return(tau)
}





