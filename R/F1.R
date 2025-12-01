########################################
# Objective function for ML estimation #
########################################
# source("R/ml_sem/lvm_ml/lvm_ml.R")


#' Compute Log-Likelihood for ML Estimation of Structural Equation Model
#'
#' This function computes the log-likelihood objective function for maximum
#' likelihood (ML) estimation of a structural equation model with latent variables.
#' It evaluates the discrepancy between the sample covariance matrix and the
#' model-implied covariance matrix.
#'
#' @param x Numeric vector containing all model parameters (loadings, correlations,
#'   path coefficients, and variance/covariance parameters).
#' @param S Sample covariance matrix of observed variables.
#' @param block_sizes Integer vector specifying the number of indicators in each block.
#' @param mode Character vector indicating the measurement mode for each block
#'   ("formative" or "reflective").
#' @param lengths_parameter Integer vector specifying the length of each parameter
#'   group in `x` (loadings, exogenous correlations, gamma, beta, endogenous correlations,
#'   variance/covariance).
#' @param which_exo_endo List containing indices and structure information for
#'   exogenous and endogenous latent variables (output from `ind_exo_endo()`).
#'
#' @return Numeric scalar representing the log-likelihood value to be minimized.
#'
#' @details
#' The log-likelihood function is computed as:
#' \deqn{\log|\Sigma(\theta)| + \text{tr}(S\Sigma(\theta)^{-1}) - \log|S| - p}
#' where:
#' - \eqn{\Sigma(\theta)} is the model-implied covariance matrix
#' - \eqn{S} is the sample covariance matrix
#' - \eqn{p} is the number of observed variables
#'
#' The function first computes the implied covariance matrix using `lvm_ml()`,
#' then evaluates the fit between sample and implied covariances.
#'
#' @seealso \code{\link{lvm_ml}} for implied covariance matrix computation
#'
#' @export
F1 <- function(x, S, block_sizes, mode, lengths_parameter, which_exo_endo){


  implied_S <- lvm_ml(x, block_sizes, mode, lengths_parameter, which_exo_endo, jac = FALSE)$SIGMA_IMPLIED

  ########################################################################
  ###################### Compute log-likelihood  #########################
  ########################################################################

  opt <- log(det(implied_S)) +
    sum(diag(S%*%solve(implied_S))) -
    log(det(S)) -
    NCOL(S)

  return(opt)

}