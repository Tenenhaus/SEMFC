


#' Compute Information Matrix for ML Estimation
#'
#' Computes the Fisher information matrix using numerical derivatives of the
#' implied covariance matrix with respect to model parameters.
#'
#' @param x Numeric vector of model parameters.
#' @param model List containing the model structure, including block sizes and
#'   other relevant information for the computation.
#'
#' @return Symmetric information matrix of dimension length(x) × length(x).
#'
#' @details
#' The information matrix is computed as:
#' \deqn{I_{ij} = \frac{1}{2} \text{tr}(\Sigma^{-1} \frac{\partial\Sigma}{\partial\theta_i} \Sigma^{-1} \frac{\partial\Sigma}{\partial\theta_j})}
#' Only the upper triangular part is computed for efficiency.
#'
#' @importFrom numDeriv jacobian
#'
#' @keywords internal
information_matrix <- function(x, model){
  JAC <- numDeriv::jacobian(lvm_ml, x = x, model = model, jac = TRUE)
  full_jac <- sapply(1:NCOL(JAC),
                        function(col){
                          ds_dt <- matrix(0, sum(model$block_sizes), sum(model$block_sizes))
                          ds_dt[upper.tri(ds_dt, diag = T)] <- JAC[, col]
                          ds_dt <- ds_dt + t(ds_dt) - diag(diag(ds_dt))
                        }, simplify = FALSE
      )

  nb_param <- length(x)
  Sinv <- solve(lvm_ml(x = x, model = model, jac = F)$SIGMA_IMPLIED)
  Sinv_full_jac <- lapply(full_jac, function(fj) as.matrix(Sinv %*% fj))
  I_ij <- function(i, j) {
    0.5 * sum(diag(Sinv_full_jac[[i]] %*% Sinv_full_jac[[j]]))
  }
  # Generate only indices of the upper triangular part
  index_upper <- which(upper.tri(matrix(0, nrow = nb_param, ncol = nb_param), diag = TRUE), arr.ind = TRUE)

  # Calculate only the elements of the upper triangular part
  I_upper <- mapply(I_ij, index_upper[, "row"], index_upper[, "col"])

  I <- matrix(0, nrow = nb_param, ncol = nb_param)
  I[upper.tri(I, diag = TRUE)] <- I_upper
  I <- I + t(I) - diag(diag(I))

  return(I)
}



#' Compute Jacobian of Constraint Functions
#'
#' Computes the transpose of the Jacobian matrix for equality constraints
#' in formative measurement models.
#'
#' @param x Numeric vector of model parameters.
#' @param S Sample covariance matrix.
#' @param model List containing the model structure, including block sizes,
#'   measurement modes, and other relevant information.
#'
#' @return Matrix H of dimension length(x) × r, where r is the number of formative blocks.
#'
#' @importFrom numDeriv jacobian
#'
#' @keywords internal
Jac_constraints <- function(x, S, model){
  # transpose of the jacobian of constraint function
  H <- t(numDeriv::jacobian(heq1, x = x, S=S, model = model))

  return(H)

}

' Compute Projection Matrix for Constrained ML Estimation
#'
#' Computes the projection matrix P for variance-covariance estimation under
#' equality constraints in structural equation models with formative blocks.
#'
#' @param x Numeric vector of model parameters.
#' @param S Sample covariance matrix.
#' @param model List containing the model structure, including:
#'   \itemize{
#'     \item \code{mode}: Character vector indicating the measurement mode for each block.
#'     \item \code{block_sizes}: Integer vector specifying the number of indicators in each block.
#'   }
#'
#' @return Projection matrix P of dimension length(x) × length(x).
#'
#' @details
#' For constrained optimization, the projection matrix is computed as:
#' \deqn{P = [M]_{1:t,1:t}^{-1}}
#' where M is the bordered information matrix:
#' \deqn{M = \begin{bmatrix} I + HH' & H \\ H' & 0 \end{bmatrix}}
#'
#'
#' @importFrom MASS ginv
#'
#' @keywords internal
P_ml <- function(x, S, model){

  mode <- model$mode

  I <- information_matrix(x, model)
  t <- nrow(I)
  r <- length(mode[mode=='formative'])
  H <- matrix(0, t, r)
  if (r>0){
    H  <- Jac_constraints(x, S, model)
  }



  M <- rbind(cbind(I+H%*%t(H), H),
             cbind(t(H), matrix(0, r, r)))

  invM <- tryCatch(
    solve(M),
    error = function(e) {
      warning("Inversion of M failed; using MASS::ginv().")
      MASS::ginv(M)
    }
  )

  P <- invM[1:t, 1:t]

  return(P)


}

#' Extract Standard Errors for Model Parameters
#'
#' Extracts and organizes standard errors for different parameter types from the
#' full standard deviation vector, including loadings, path coefficients, regression
#' coefficients, effects, and residual variances.
#'
#' @param SD Numeric vector of standard deviations for all model parameters.
#' @param mode Character vector indicating the measurement mode for each block
#'   ("reflective" or "formative").
#' @param lengths_parameter Integer vector specifying the length of each parameter group.
#' @param block_sizes Integer vector specifying the number of indicators in each block.
#' @param vcov_effect List containing variance-covariance matrices for effects:
#'   \itemize{
#'     \item \code{vcov_exo_total}: Variance-covariance matrix for total effects on exogenous variables.
#'     \item \code{vcov_endo_total}: Variance-covariance matrix for total effects on endogenous variables.
#'     \item \code{vcov_exo_indirect}: Variance-covariance matrix for indirect effects on exogenous variables.
#'     \item \code{vcov_endo_indirect}: Variance-covariance matrix for indirect effects on endogenous variables.
#'   }
#'
#' @return List containing standard errors for each parameter type:
#'   \item{sd_lambda}{Standard errors for loadings.}
#'   \item{sd_gamma}{Standard errors for path coefficients between latent variables.}
#'   \item{sd_beta}{Standard errors for regression coefficients.}
#'   \item{sd_total_effects}{Standard errors for total effects.}
#'   \item{sd_indirect_effects}{Standard errors for indirect effects.}
#'   \item{sd_residual_variance}{Standard errors for residual variances (reflective blocks only).}
#'
#' @keywords internal
get_se_series <- function(SD, mode, lengths_parameter, block_sizes, vcov_effect){
  start_indices_in_x <- cumsum(c(1, head(lengths_parameter, -1)))
  lambda_start_index <- start_indices_in_x[1]
  lambda_end_index <- lambda_start_index + lengths_parameter[1] - 1
  gamma_start_index <- start_indices_in_x[3]
  gamma_end_index <- gamma_start_index +  lengths_parameter[3] - 1
  beta_start_index <- start_indices_in_x[4]
  beta_end_index <- beta_start_index + lengths_parameter[4] - 1

  sd_lambda <- SD[lambda_start_index: lambda_end_index]
  sd_gamma <- SD[gamma_start_index: gamma_end_index]
  sd_beta <- SD[beta_start_index: beta_end_index]
  sd_total_effects <- sqrt(c(diag(vcov_effect$vcov_exo_total), diag(vcov_effect$vcov_endo_total)))
  sd_indirect_effects <- sqrt(c(diag(vcov_effect$vcov_exo_indirect), diag(vcov_effect$vcov_endo_indirect)))


  BDIAG <- get_bdiag(SD,
                     mode = mode,
                     block_sizes = block_sizes,
                     initial_start_index_cov = start_indices_in_x[6])
  sd_residual_variance <- unlist(lapply(BDIAG[mode == 'reflective'], diag))

  return(list(
    sd_lambda = sd_lambda,
    sd_gamma = sd_gamma,
    sd_beta = sd_beta,
    sd_total_effects = sd_total_effects,
    sd_indirect_effects = sd_indirect_effects,
    sd_residual_variance = sd_residual_variance
  ))


}


#' Format ML Inference Results
#'
#' Organizes parameter estimates, standard errors, z-scores and p-values into
#' structured data frames for loadings, path coefficients, regression coefficients,
#' effects, and residual variances.
#'
#' @param fit List containing model fit results from `lvm_ml()`. Includes parameter estimates
#'   such as loadings, standardized loadings, gamma coefficients, beta coefficients,
#'   residual variances, and effects.
#' @param model List containing the model structure, including:
#'   \itemize{
#'     \item \code{mode}: Character vector indicating the measurement mode for each block.
#'     \item \code{lengths_theta}: Integer vector specifying the length of each parameter group.
#'     \item \code{block_sizes}: Integer vector specifying the number of indicators in each block.
#'   }
#' @param VCOV Variance-covariance matrix of parameter estimates.
#' @param vcov_effect List containing variance-covariance matrices for effects:
#'   \itemize{
#'     \item \code{vcov_exo_total}: Variance-covariance matrix for total effects on exogenous variables.
#'     \item \code{vcov_endo_total}: Variance-covariance matrix for total effects on endogenous variables.
#'     \item \code{vcov_exo_indirect}: Variance-covariance matrix for indirect effects on exogenous variables.
#'     \item \code{vcov_endo_indirect}: Variance-covariance matrix for indirect effects on endogenous variables.
#'   }
#'
#' @return List of data frames with columns: lhs, op, rhs, est, se, z, ci.lower, ci.upper, pvalue, std.all:
#'   \item{lambda}{Loadings estimates and inference statistics.}
#'   \item{gamma}{Path coefficients between latent variables.}
#'   \item{beta}{Regression coefficients.}
#'   \item{total_effects}{Total effects.}
#'   \item{indirect_effects}{Indirect effects.}
#'   \item{residual_variance}{Residual variances (reflective blocks only).}
#'
#' @keywords internal
formatting_ml_infer <- function(fit, model, VCOV, vcov_effect ){


  mode <- model$mode
  lengths_parameter <- model$lengths_theta
  block_sizes <- model$block_sizes
  SD <- sqrt(diag(VCOV))


  se <- get_se_series(SD, mode, lengths_parameter, block_sizes, vcov_effect)
  table <- formatting_estimate(fit, se)

  return(table)


}


#' Statistical Inference for ML-SEM
#'
#' Performs statistical inference for maximum likelihood structural equation models,
#' computing variance-covariance matrix and organizing results into tables.
#'
#' @param x Numeric vector of optimal parameter values.
#' @param S Sample covariance matrix.
#' @param model List containing the model structure, including:
#'   \itemize{
#'     \item \code{mode}: Character vector indicating the measurement mode for each block.
#'     \item \code{block_sizes}: Integer vector specifying the number of indicators in each block.
#'     \item \code{lengths_theta}: Integer vector specifying the length of each parameter group.
#'   }
#' @param N Integer sample size.
#' @param fit List containing model fit results from `lvm_ml()`, including parameter estimates
#'   such as loadings, path coefficients, and residual variances.
#'
#' @return List containing:
#'   \item{estimate}{List of data frames with parameter estimates and inference statistics.}
#'   \item{VCOV}{Variance-covariance matrix of parameter estimates.}
#'   \item{SD}{Vector of standard errors for all parameters.}
#'
#' @details
#' The variance-covariance matrix is computed as VCOV = P/N, where P is the
#' projection matrix and N is the sample size.
#'
#' @keywords internal
mlSEM_infer <- function(x, S, model, N, fit){

  P_ml <- P_ml(x, S, model)
  VCOV <- P_ml/N
  # SD <- sqrt(diag(VCOV))

  vcov_effect <- effect_infer(fit$beta, fit$gamma, model$lengths_theta, VCOV)

  table <- formatting_ml_infer(fit, model, VCOV, vcov_effect)

  out <- list(
    estimate = table,
    VCOV = VCOV,
    vcov_effect = vcov_effect

  )

  return(out)


}


#' Compute Z-Score and P-Value for Hypothesis Test
#'
#' Computes the z-score and p-value for testing linear hypotheses H0: L*theta = 0
#' on model parameters.
#'
#' @param x Numeric vector of optimal parameter values.
#' @param S Sample covariance matrix.
#' @param X Data matrix (currently unused in function body).
#' @param model List containing the model structure, including:
#'   \itemize{
#'     \item \code{mode}: Character vector indicating the measurement mode for each block.
#'     \item \code{block_sizes}: Integer vector specifying the number of indicators in each block.
#'   }
#' @param L Contrast matrix defining the linear hypothesis.
#'
#' @return P-value for the two-tailed test.
#'
#' @details
#' The z-score is computed as:
#' \deqn{z = \frac{\sqrt{N} L\theta}{\sqrt{LPL'}}}
#' The p-value is two-tailed using the standard normal distribution.
#'
#' @keywords internal
z_H0 <- function(x, S, X, model, L){
  N <- nrow(X)
  P <- P_ml(x, S, model)
  z <- sqrt(N)*(L %*% x)/sqrt(L%*%P%*%t(L))
  pval <- 2*pnorm(z, lower.tail = F)

  return(pval)

}















