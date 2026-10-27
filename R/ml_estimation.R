

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
#' @param constraint_tol Numeric tolerance for checking the maximum
#'   absolute equality-constraint violation (default: 1e-6).
#'
#' @return A list containing:
#' \describe{
#'   \item{pars}{Numeric vector of estimated parameters.}
#'   \item{diagnostics}{List containing the solver name, convergence
#'     status, convergence code, message, objective value, maximum
#'     constraint violation, constraint satisfaction indicator,
#'     and overall success indicator.}
#'   \item{raw_result}{Original result returned by the optimizer
#'     (`nloptr` or `nlminb`).}
#' }
#'
#' @details
#' The function minimizes the maximum-likelihood fitting function
#' `F1()` using analytical derivatives. The optimization method
#' depends on the presence of formative blocks:
#'
#' - Models with formative blocks: Constrained optimization is
#'   performed using the SLSQP algorithm (`NLOPT_LD_SLSQP`) from
#'   the `nloptr` package. Equality constraints defined by `heq()`
#'   enforce the identification conditions of the composite blocks.
#'   Analytical gradients of the objective function and the
#'   constraints are provided through `grad_F1()` and `grad_heq()`.
#'
#' - Models without formative blocks: Unconstrained optimization
#'   is performed using `nlminb()`, with analytical gradient and
#'   Hessian functions provided by `grad_F1()` and `hess_F1()`.
#'
#' The optimization tolerance is controlled by `tol`. For SLSQP,
#' the maximum number of function evaluations is set to 500.
#'
#' Convergence and constraint satisfaction are assessed separately.
#' A solution is considered successful when the optimizer reports
#' convergence, the objective value is finite, and the maximum
#' absolute constraint violation does not exceed `constraint_tol`.
#' The function returns the estimated parameters together with
#' optimization diagnostics and the raw solver output.
#'
#' @importFrom nloptr nloptr
#' @importFrom stats nlminb
#'
#' @keywords internal
mlSEM <- function(init, S, model, tol = 1e-8,
                  constraint_tol = 1e-6) {

  mode <- model$mode
  r <- sum(mode == "formative")

  solver <- if (r > 0L) "SLSQP" else "nlminb"
  f <- function (x) { return(F1(x, S, model)) }
  grad <- function(x){ return(grad_F1(x, S, model)) }
  hess <- function(x){ return(hess_F1(x, S, model)) }
  if (r > 0L) {
    h_eq <- function(x){ return(heq(x, S, model)) }
    grad_h_eq <- function(x){ return(grad_heq(x, S, model)) }
  }
  result <- tryCatch({

    if (r > 0L) {
      nloptr::nloptr(
        x0 = init,
        eval_f = f,
        eval_grad_f = grad,
        eval_g_eq = h_eq,
        eval_jac_g_eq = grad_h_eq,
        opts = list(
          algorithm = "NLOPT_LD_SLSQP",
          xtol_rel = tol,
          tol_constraints_eq = rep(constraint_tol, r),
          maxeval = 500L
        )
      )

    } else {

      stats::nlminb(
        start = init,
        objective = f,
        gradient = grad,
        hessian = hess,
        control = list(x.tol = tol)
      )

    }

  }, error = function(e) e)

  if (inherits(result, "error")) {

    diagnostics <- list(
      solver = solver,
      converged = FALSE,
      convergence_code = NA_integer_,
      message = conditionMessage(result),
      objective = NA_real_,
      max_constraint_violation = if (r > 0L) NA_real_ else 0,
      constraints_satisfied = (r == 0L),
      success = FALSE
    )

    return(list(
      pars = init,
      diagnostics = diagnostics,
      raw_result = NULL
    ))
  }
  if (r > 0L) {
    x <- result$solution
    objective <- result$objective
    convergence_code <- result$status
    message <- result$message
    converged <- isTRUE(convergence_code %in% 1:4)

  } else {
    x <- result$par
    objective <- result$objective
    convergence_code <- result$convergence
    message <- result$message
    converged <- isTRUE(convergence_code == 0)
  }
  if (r > 0L) {
    max_violation <- max(abs(h_eq(x)))
  } else {
    max_violation <- 0
  }

  constraints_satisfied <-
    is.finite(max_violation) &&
    max_violation <= constraint_tol

  diagnostics <- list(
    solver = solver,
    converged = converged,
    convergence_code = convergence_code,
    message = message,
    objective = objective,
    max_constraint_violation = max_violation,
    constraints_satisfied = constraints_satisfied,
    success = converged &&
              constraints_satisfied &&
              isTRUE(is.finite(objective))
  )
  return(list(
    pars = x,
    diagnostics = diagnostics,
    raw_result = result
  ))
}


#' Compute Log-Likelihood for ML Estimation of Structural Equation Model
#'
#' This function computes the log-likelihood objective function for maximum
#' likelihood estimation of a structural equation model with latent variables.
#' It evaluates the discrepancy between the sample covariance matrix and the
#' model-implied covariance matrix.
#'
#' @param x Numeric vector containing all model parameters (loadings, correlations,
#'   path coefficients, and variance/covariance parameters).
#' @param S Sample covariance matrix of observed variables.
#' @param model A list containing model specifications with the following elements:
#'   \describe{
#'     \item{block_sizes}{Integer vector specifying the number of indicators in each block.}
#'     \item{mode}{Character vector indicating the measurement mode for each block
#'       ("formative" or "reflective").}
#'     \item{lengths_parameter}{Integer vector specifying the length of each parameter
#'       group in `x` (loadings, exogenous correlations, gamma, beta, endogenous correlations,
#'       covariance of composite blocks and residual variances).}
#'     \item{which_exo_endo}{List containing indices and structure information for
#'       exogenous and endogenous latent variables (output from `ind_exo_endo()`).}
#'    \item{dag}{Logical indicating whether the structural model is recursive (FALSE)
#'       or non-recursive (TRUE).}
#'   }
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
#' @seealso \code{\link{lvm_ml}} for implied covariance matrix computation,
#'   \code{\link{ind_exo_endo}} for structure information generation
#'
#'
#' @examples
#' \dontrun{
#' # Example with 6 blocks: 4 formative exogenous, 2 reflective endogenous
#' set.seed(123)
#' n <- 300
#'
#'
#' # Generate sample data (6 blocks with 3 indicators each)
#' X <- matrix(rnorm(n * 18), ncol = 18)
#' S <- cov(X)
#'
#' # Model specification
#' block_sizes <- c(3, 3, 3, 3, 3, 3)
#' mode <- c("formative", "formative", "formative", "formative",
#'           "reflective", "reflective")
#'
#' # Structure information (hard-coded for this example)
#' which_exo_endo <- list(
#'   Hi = list(
#'     c(LV1 = 1, LV2 = 2),
#'     c(LV3 = 3, LV4 = 4)
#'   ),
#'   Ji = list(
#'     c(LV6 = 6),
#'     c(LV5 = 5)
#'   ),
#'   ind_endo = c(LV5 = 5, LV6 = 6),
#'   ind_exo = c(LV1 = 1, LV2 = 2, LV3 = 3, LV4 = 4)
#' )
#'
#' dag <- FALSE
#'
#' varnames <- list(
#'  LV1 = c("X11", "X12", "X13"),
#'  LV2 = c("X21", "X22", "X23"),
#'  LV3 = c("X31", "X32", "X33"),
#'  LV4 = c("X41", "X42", "X43"),
#'  LV5 = c("X51", "X52", "X53"),
#'  LV6 = c("X61", "X62", "X63")
#')
#'
#'
#' # Parameter vector structure:
#' # - 18 loadings (3 per  block)
#' # - 6 exogenous correlations
#' # - 4 non zero gamma coefficients
#' # - 2 non zero beta coefficient
#' # - 1 endogenous correlations
#' # - 30 covariances and residual variances
#' lengths_theta <- c(18, 6, 4, 2, 1, 30)
#'
#'
#' model <- list(
#'  block_sizes = block_sizes,
#' mode = mode,
#' varnames = varnames,
#' lengths_theta = lengths_theta,
#' which_exo_endo = which_exo_endo,
#' dag = dag
#' )
#'
#'
#'
#' # Initialize parameters
#' x  <- rnorm(61)
#' # Compute log-likelihood
#' f <- F1(x, S, model)
#' print(f)
#'}
#' @keywords internal
F1 <- function(x, S, model){


  implied_S <- lvm_ml(x, model, jac = FALSE)$sigma_implied

  # check for positive definiteness

  eigvals <- eigen(implied_S, symmetric = TRUE, only.values = TRUE)$values
  if (any(eigvals <= 0.0001)){
    return(1e+5)
  }

  opt <- log(det(implied_S)) +
    sum(diag(S%*%solve(implied_S))) -
    log(det(S)) -
    NCOL(S)

  return(opt)

}



#' @title Equality Constraints for Loglikelihood Minimization (H1)
#'
#' @description
#' Generalized version for multiple blocks with formative/reflective modes.
#'
#' @param x Numeric vector of parameters.
#' @param S Sample covariance matrix.
#' @param model A list containing model specifications with the following elements:
#'   \describe{
#'     \item{block_sizes}{Integer vector specifying the number of indicators in each block.}
#'     \item{mode}{Character vector indicating the measurement mode for each block
#'       ("formative" or "reflective").}
#'   }
#'
#' @return Vector of equality constraints.
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' n <- 300
#' # Generate sample data (6 blocks with 3 indicators each)
#' X <- matrix(rnorm(n * 18), ncol = 18)
#' S <- cov(X)
#' # Model specification
#' block_sizes <- c(3, 3, 3, 3, 3, 3)
#' mode <- c("formative", "formative", "formative", "formative",
#'           "reflective", "reflective")
#'
#'
#' model <- list(
#'  block_sizes = block_sizes,
#' mode = mode
#' )
#'
#' x  <- rnorm(61)
#'
#' heq1(x, S, model)
#'}
#' @keywords internal
heq1 <- function (x, S, model){

  block_sizes <- model$block_sizes
  mode <- model$mode

  loadings <- get_loadings(x, block_sizes)



 # list of lengths: diagonal covariances for reflective blocks,
 # or lower-triangular covariance values for formative blocks
  lengths_values_cov <- block_sizes
  lengths_values_cov[mode == "formative"] <- (block_sizes[mode == "formative"]^2 + block_sizes[mode == "formative"]) / 2
  # number of parameters for covariance
  total_cov_parameter <- sum(lengths_values_cov)
    # the coefficient are stocked at the end of x
  initial_start_index_cov <- length(x) - total_cov_parameter +1
  end_endex_cov <- initial_start_index_cov + total_cov_parameter - 1

  # part of the vector corresponding to covariance blocks
  extracted_parameters_cov <- x[initial_start_index_cov:end_endex_cov]
  # list of parameters corresponding to each covariance bloc
  list_cov <- split(extracted_parameters_cov,
                    rep(seq_along(lengths_values_cov), lengths_values_cov))

  # list of formative covariance matrices
  S_composite <- lapply(list_cov[mode == "formative"], build_formative_S_diag)
  h <- as.vector(mapply(function(S_composite_i, loadings_i) {
    t(loadings_i) %*% solve(S_composite_i) %*% loadings_i - 1
  }, S_composite, loadings[mode == "formative"], SIMPLIFY = TRUE))

  return(h)

}




heq <- function (x, S, model){

  block_sizes <- model$block_sizes
  mode <- model$mode
  r <- sum(mode == "formative")

  loadings <- get_loadings(x, block_sizes)
  loadings_formative <- loadings[mode == "formative"]
  Lambda_formative <- Matrix::bdiag(loadings_formative)



 # list of lengths: diagonal covariances for reflective blocks,
 # or lower-triangular covariance values for formative blocks
  lengths_values_cov <- block_sizes
  lengths_values_cov[mode == "formative"] <- (block_sizes[mode == "formative"]^2 + block_sizes[mode == "formative"]) / 2
  # number of parameters for covariance
  total_cov_parameter <- sum(lengths_values_cov)
    # the coefficient are stocked at the end of x
  initial_start_index_cov <- length(x) - total_cov_parameter +1
  end_endex_cov <- initial_start_index_cov + total_cov_parameter - 1

  # part of the vector corresponding to covariance blocks
  extracted_parameters_cov <- x[initial_start_index_cov:end_endex_cov]
  # list of parameters corresponding to each covariance bloc
  list_cov <- split(extracted_parameters_cov,
                    rep(seq_along(lengths_values_cov), lengths_values_cov))

  # list of formative covariance matrices
  S_composite <- lapply(list_cov[mode == "formative"], build_formative_S_diag)


  inv_Sigma_formative <- bdiag(lapply(S_composite, solve))

  h <- diag(t(Lambda_formative) %*% inv_Sigma_formative %*% Lambda_formative) - rep(1, r)

  return(h)

}


grad_heq <- function(x, S, model) {
  est <- lvm_ml(x, model, jac = FALSE)
  lambda <- est$lambda
  S_composites <- est$S_composites


  J_H <- compute_gradient_constraint(lambda, S_composites, model$block_sizes, model$lengths_theta, model$mode)

  return(as.matrix(J_H))
}

