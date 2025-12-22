

# library(Rsolnp)
#
# source('R/ml_sem/F1.R')
# source('R/ml_sem/h_constraints.R')


#' Maximum Likelihood Estimation of Structural Equation Model
#'
#' Performs ML estimation of a structural equation model using constrained
#' optimization via the solnp algorithm. Handles both formative and reflective
#' measurement models with automatic constraint application for formative blocks.
#'
#' @param init Numeric vector of initial parameter values for optimization.
#' @param block_sizes Integer vector specifying the number of indicators in each block.
#' @param mode Character vector indicating the measurement mode for each block
#'   ("formative" or "reflective").
#' @param S Sample covariance matrix of observed variables.
#' @param lengths_parameter Integer vector specifying the length of each parameter
#'   group (loadings, exogenous correlations, gamma, beta, endogenous correlations,
#'   variance/covariance).
#' @param which_exo_endo List containing indices and structure information for
#'   exogenous and endogenous latent variables (output from `ind_exo_endo()`).
#'
#' @return Object of class "solnp" containing optimization results:
#'   \item{pars}{Optimal parameter values.}
#'   \item{convergence}{Convergence code (0 = success).}
#'   \item{values}{Objective function values at each iteration.}
#'   \item{hessian}{Hessian matrix at optimum.}
#'
#' @details
#' The function uses `solnp()` from the Rsolnp package to minimize the
#' log-likelihood function `F1()`. When formative blocks are present (r > 0),
#' equality constraints (`heq1()`) are automatically applied to ensure proper
#' identification of the model. The optimization uses a tolerance of 1e-4
#' with trace output disabled.
#'
#' @importFrom Rsolnp solnp
#'
#' @export
mlSEM <- function (init, S, model){

  mode <- model$mode
  # number of formative blocks
  r <- sum(mode == "formative")

  if(r !=0){

  result <- solnp(pars = init,
                  fun=F1, eqfun=heq1,
                  eqB = rep(0,r), S = S, model = model,
                  control = list(trace = 0, tol = 1e-4))
  }
  else{
    result <- solnp(pars = init,
                    fun=F1, S = S, model = model,
                    control = list(trace = 0, tol = 1e-4))

  }

  return(result)



}