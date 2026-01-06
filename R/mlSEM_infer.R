


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
information_matrix <- function(x, model, data){
  JAC <- numDeriv::jacobian(lvm_ml, x = x, model = model, data = data, jac = TRUE)
  full_jac <- sapply(1:NCOL(JAC),
                        function(col){
                          ds_dt <- matrix(0, sum(model$block_sizes), sum(model$block_sizes))
                          ds_dt[upper.tri(ds_dt, diag = T)] <- JAC[, col]
                          ds_dt <- ds_dt + t(ds_dt) - diag(diag(ds_dt))
                        }, simplify = FALSE
      )

  nb_param <- length(x)
  Sinv <- solve(lvm_ml(x = x, model = model, data =  data)$SIGMA_IMPLIED)
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


#' Format ML Inference Results
#'
#' Organizes parameter estimates, standard errors, z-scores and p-values into
#' structured data frames for loadings, path coefficients and residual variances.
#'
#' @param fit List containing model fit results from `lvm_ml()`. Includes parameter estimates
#'   such as loadings, standardized loadings, gamma coefficients, beta coefficients, and residual variances.
#' @param SD Numeric vector of standard errors for all parameters.
#' @param model List containing the model structure, including:
#'   \itemize{
#'     \item \code{mode}: Character vector indicating the measurement mode for each block.
#'     \item \code{lengths_theta}: Integer vector specifying the length of each parameter group.
#'     \item \code{block_sizes}: Integer vector specifying the number of indicators in each block.
#'   }
#'
#' @return List containing:
#'   \item{lambda}{Data frame with loadings estimates, std errors, z-scores and p-values.}
#'   \item{std_lambda}{Data frame with standardized loadings estimates, std errors, z-scores and p-values.}
#'   \item{gamma}{Data frame with gamma coefficients estimates, std errors, z-scores and p-values.}
#'   \item{beta}{Data frame with beta coefficients estimates, std errors, z-scores and p-values.}
#'   \item{residual_variance}{Data frame with residual variances estimates, std errors, z-scores and p-values.}
#'
#' @details
#' Z-scores are computed as estimate/std error. P-values are two-tailed using
#' the standard normal distribution.
#'
#' @keywords internal
formatting_ml_infer <- function(fit, SD, model){


  mode <- model$mode
  lengths_parameter <- model$lengths_theta
  block_sizes <- model$block_sizes


  lambda <- unlist(fit$lambda)
  std_lambda <- unlist(fit$std_lambda)
  gamma <- fit$gamma[fit$gamma!=0]
  beta <- fit$beta[fit$beta!=0]
  residual_variance <- unlist(unname(fit$residual_variance))



  start_indices_in_x <- cumsum(c(1, head(lengths_parameter, -1)))
  lambda_start_index <- start_indices_in_x[1]
  lambda_end_index <- lambda_start_index + length(lambda) - 1
  gamma_start_index <- start_indices_in_x[3]
  gamma_end_index <- gamma_start_index + length(gamma) - 1
  beta_start_index <- start_indices_in_x[4]
  beta_end_index <- beta_start_index + length(beta) - 1

  sd_lambda <- SD[lambda_start_index: lambda_end_index]
  sd_std_lambda <- sd_lambda*std_lambda/lambda
  sd_gamma <- SD[gamma_start_index: gamma_end_index]
  sd_beta <- SD[beta_start_index: beta_end_index]

  BDIAG <- get_bdiag(SD,
                     mode = mode,
                     block_sizes = block_sizes,
                     initial_start_index_cov = start_indices_in_x[6])
  sd_residual_variance <- unlist(lapply(BDIAG[mode == 'reflective'], diag))



  z_lambda <- lambda/sd_lambda
  z_gamma<- gamma/sd_gamma
  z_beta <- beta/sd_beta
  z_residual_variance <- residual_variance/sd_residual_variance


  table_lambda <- data.frame(Estimate = lambda,
                             std = sd_lambda,
                             z_score = z_lambda,
                             ci_lower = lambda - 1.96 * sd_lambda,
                             ci_upper = lambda + 1.96 * sd_lambda,
                             pval = unlist(lapply(z_lambda, function (z) 2*pnorm(abs(z), lower.tail = FALSE)))

  )
  rownames(table_lambda) <- gsub("\\.", "~", rownames(table_lambda))


  table_std_lambda <- data.frame(Estimate = std_lambda,
                           std = sd_std_lambda,
                           z_score = z_lambda,
                           ci_lower = std_lambda - 1.96 * sd_std_lambda,
                           ci_upper = std_lambda + 1.96 * sd_std_lambda,
                           pval = unlist(lapply(z_lambda, function (z) 2*pnorm(abs(z), lower.tail = FALSE)))

  )
  rownames(table_std_lambda) <- rownames(table_lambda)







  table_gamma <- data.frame(Estimate = gamma,
                            std = sd_gamma,
                            z_score = z_gamma,
                            ci_lower = gamma - 1.96 * sd_gamma,
                            ci_upper = gamma + 1.96 * sd_gamma,

                            pval = unlist((lapply(z_gamma, function (z) 2*pnorm(abs(z), lower.tail = FALSE))))
  )

  rownames(table_gamma) <- sapply(1:NROW(table_gamma),
                          function(b)
                            paste(rownames(fit$gamma)[which(fit$gamma!=0, arr.ind = TRUE)[b, 1]],
                                  colnames(fit$gamma)[which(fit$gamma!=0, arr.ind = TRUE)[b, 2]],
                                  sep = "~")
  )

  table_beta <- data.frame()

  if (length(beta) != 0){

  table_beta <- data.frame(Estimate = beta,
                           std = sd_beta,
                           z_score = z_beta,
                           ci_lower = beta - 1.96 * sd_beta,
                           ci_upper = beta + 1.96 * sd_beta,
                           pval = unlist(lapply(z_beta, function (z) 2*pnorm(abs(z), lower.tail = FALSE)))
  )
  rownames(table_beta) <- sapply(1:NROW(table_beta),
       function(b)
         paste(colnames(fit$beta)[which(fit$beta!=0, arr.ind = TRUE)[b, ]],
               collapse = "~")
       )
  }

  table_residual_variance <- data.frame(Estimate = residual_variance,
                           std = sd_residual_variance,
                           z_score = z_residual_variance,
                           ci_lower = residual_variance - 1.96 * sd_residual_variance,
                           ci_upper = residual_variance + 1.96 * sd_residual_variance,
                           pval = unlist(lapply(z_residual_variance, function (z) 2*pnorm(abs(z), lower.tail = FALSE)))
  )



  out <- list(
    lambda = table_lambda,
    std_lambda = table_std_lambda,
    gamma = table_gamma,
    beta = table_beta,
    residual_variance = table_residual_variance

  )

  return(out)


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
  SD <- sqrt(diag(VCOV))

  table <- formatting_ml_infer(fit, SD, model)

  out <- list(
    estimate = table,
    VCOV = VCOV,
    SD = SD
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















