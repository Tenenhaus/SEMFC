

#' Global Fit Assessment for svdSEM Models
#'
#' This function computes the global fit test statistic for a svdSEM model
#' and optionally performs a bootstrap test to assess the significance of the fit.
#' The method is based on the Yuan & Hayashi (2003) approach for data transformation.
#'
#' @param fit A svdSEM fit object containing the fitted model information
#' @param B Number of bootstrap replications. If NULL, only the T_LS statistic is returned
#' @param bias Logical indicating whether a biased covariance estimator should be used (default: FALSE)
#'
#' @return If B is NULL, returns the T_LS statistic. Otherwise, returns a list containing:
#'   \item{T_LS}{The observed test statistic}
#'   \item{Tb_LS}{Vector of bootstrap test statistics}
#'   \item{pval}{The p-value of the fit test}
#'   \item{improper}{Number of improper solutions (negative eigenvalues) encountered}
#'
#' @details The function first transforms the data according to Yuan & Hayashi (2003),
#' then performs bootstrap to estimate the empirical distribution of the test statistic.
#' Solutions with negative eigenvalues are considered improper and excluded.
#'
#' @examples
#' \dontrun{
#' gof_results <- svdSEM_gof(fit, B = 100, bias = FALSE)
#' print(gof_results$pval)
#' }
#'
#'
#' @keywords internal

svdSEM_gof <- function(fit, B = 100, bias = FALSE){
  
  if(is.null(B)){
    return(T_LS = fit$T_LS)
  }else{
    # Tranforms the data sets in the way proposed 
    # by Yuan & Hayashi (2003)
    df <- Reduce("cbind", fit$blocks)
    Z0 <- scaleDataSet(df, fit$sigma_implied)
    Z0 <- lapply(split(data.frame(t(Z0)),
                      as.factor(rep(seq_along(fit$blocks),
                                    sapply(fit$blocks, NCOL)))),
                t)
    
    beta = fit$beta
    gamma = fit$gamma
    Tb_LS = rep(NA, B)
    
    Tb_LS <- pbapply::pbsapply(1:B,
                      function(b){
                        ind <- sample(NROW(df), replace = TRUE)
                        Zb <- lapply(Z0, function(x) x[ind, ])
                        S <- cov2(Reduce("cbind", Zb), bias = bias)
                        fit_b <- svdSEM(Zb,
                                       C = fit$C,
                                       scale = fit$scale, 
                                       mode = fit$mode, 
                                       bias = fit$bias)
                        if (!any(eigen(fit_b$p_tilde)$values<=0)){
                          return(NA)
                        }else{
                          return(d_LS(S, fit_b$sigma_implied))
                        }
                      }
                      )
    
    
    pval <- mean(fit$T_LS <= Tb_LS, na.rm = TRUE)
    return(list(T_LS = fit$T_LS,
                Tb_LS = Tb_LS,
                pval = pval, 
                improper = sum(is.na(Tb_LS))))
  }
}





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
#' @keywords internal

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

  

