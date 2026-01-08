

#' Maximum Likelihood Estimation of Structural Equation Model
#'
#' Performs ML estimation of a structural equation model using constrained
#' optimization via the solnp algorithm. Handles both formative and reflective
#' measurement models with automatic constraint application for formative blocks.
#'
#' @param init Numeric vector of initial parameter values for optimization.
#' @param S Sample covariance matrix of observed variables.
#' @param model A list containing model specifications with the following elements:
#'   \describe{
#'     \item{mode}{Character vector indicating the measurement mode for each block}
#'   }
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
#' identification of the model. The optimization uses a tolerance of 1e-8
#' with trace output disabled.
#'
#' @importFrom Rsolnp solnp
#'
#' @export
mlSEM <- function (init, S, model, tol = 1e-08){

  mode <- model$mode
  # number of formative blocks
  r <- sum(mode == "formative")

  if(r !=0){

  result <- solnp(pars = init,
                  fun=F1, eqfun=heq1,
                  eqB = rep(0,r), S = S, model = model,
                  control = list(trace = 0, tol = tol))
  }
  else{
    result <- solnp(pars = init,
                    fun=F1, S = S, model = model,
                    control = list(trace = 0, tol = tol))

  }

  return(result)



}