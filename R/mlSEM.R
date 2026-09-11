

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
#' @param tol Numeric value specifying the optimization tolerance (default: 1e-8).
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
#' @keywords internal
mlSEM <- function (init, S, model, tol = 1e-08){

  mode <- model$mode
  # number of formative blocks
  r <- sum(mode == "formative")

  f <- function (x) { return(F1(x, S, model)) }
  grad <- function(x){ return(grad_F1(x, S, model)) }
  h_eq <- function(x){ return(heq(x, S, model)) }
  grad_h_eq <- function(x){ return(grad_heq(x, S, model)) }
  hess <- function(x){ return(hess_F1(x, S, model)) }

  if(r !=0){
    result <- nloptr(
      x0=init,
      eval_f=f,
      eval_grad_f=grad,
      eval_g_eq = h_eq,
      eval_jac_g_eq = grad_h_eq,
      opts = list("algorithm"="NLOPT_LD_SLSQP",
                  'xtol_rel' = tol,
                  'maxeval' = 500)
    )
    x <- result$solution
  }
  else{
    result <- nlminb(
      start = init,
      objective = f,
      gradient = grad,
      hessian = hess,
      lower = -Inf,
      upper = Inf,
      control = list(x.tol = tol)
    )
    x <- result$par

  }

  return(x)



}