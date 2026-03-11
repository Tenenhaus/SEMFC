
#' Detect Improper Solutions in SEM Estimation
#'
#' Checks for nine types of inadmissible or improper solutions in structural
#' equation models, including out-of-range estimates and non-positive definite
#' matrices.
#'
#' @param fit A list of estimates, including:
#' \itemize{
#'   \item \code{reliability_coef}: Reliability coefficients for each block
#'   \item \code{p_implied}: Implied correlation matrix of latent variables
#'   \item \code{p_tilde}: Corrected correlation matrix of latent variables
#'   \item \code{residual_variance}: List of residual variances
#'   \item \code{std_lambda}: Standardized loadings
#'   \item \code{sigma_implied}: Implied covariance matrix of observed variables
#'   \item \code{r2}: R-squared values for endogenous latent variables
#'   \item \code{psi}: Residual covariance matrix of latent variables
#' }
#'
#' @return Named logical vector of length 8 (or 9 for SVDSEM) indicating presence of each type
#'   of improper solution:
#' \describe{
#'   \item{reliability_coef}{\code{TRUE} if any reliability coefficient is outside (0, 1).}
#'   \item{rho_jh}{\code{TRUE} if any correlation in p_implied is outside (-1, 1).}
#'   \item{p_implied}{\code{TRUE} if p_implied has negative eigenvalues (not positive definite).}
#'   \item{theta_jh}{\code{TRUE} if any residual variance is negative.}
#'   \item{std_lambda}{\code{TRUE} if any standardized loading is outside (-1, 1).}
#'   \item{sigma_implied}{\code{TRUE} if sigma_implied has negative eigenvalues (not positive definite).}
#'   \item{r2}{\code{TRUE} if any R-squared is outside (0, 1).}
#'   \item{psi}{\code{TRUE} if psi has negative eigenvalues (not positive definite).}
#'   \item{p_tilde}{\code{TRUE} if p_tilde has negative eigenvalues (only for SVDSEM).}
#' }
#'
#' @details
#' A value of \code{TRUE} indicates an improper solution for that criterion. An
#' admissible solution should have all values set to \code{FALSE}. Common causes
#' of improper solutions include poor model specification, insufficient sample
#' size, or convergence to a boundary solution.
#'
#' @keywords internal
improper <- function(fit){

  
  improper_sol <- rep(FALSE, 8)
  
  # Improper solution 1 : number of out of range reliability coefficient
  # Expected value : 0 <=rho_jk <=1

  improper_sol[1] <- any(fit$reliability_coef < 0 |
                        fit$reliability_coef > 1)
  
  
  # Improper solution 2 : number of out of range rho_jk 
  # Expected value : -1 <=rho_jk <=1
  
  improper_sol[2] <- any(abs(fit$p_implied[upper.tri(fit$p_implied)])>1)
  

  
  # Improper solution 3 : is p_implied positive definite
  # Expected value = TRUE
  
  improper_sol[3] <- any(eigen(fit$p_implied)$values<0)
  
  # Improper solution 4 : is variance of the residual positive
  # Expected value = TRUE
  
  improper_sol[4] <- any(unlist(fit$residual_variance)<0)
  
  # Improper solution 5 : number of out of range cor(y_jh, eta_j)
  # Expected value : -1 <= cor(y_jh, eta_j) <=1
  
  improper_sol[5] <- any(abs(unlist(fit$std_lambda))>1)
  
  # Improper solution 6 : is sigma_implied positive definite
  # Expected value = TRUE
  
  improper_sol[6] <- any(eigen(fit$sigma_implied)$values<0)
  
  # Improper solution 7 : number of out of range R^2_i
  # Expected value : 0 <= R^2_i <=1
  
  improper_sol[7] <- any(fit$r2< 0 | fit$r2>1)
  
  # Improper solution 8 : is psi positive definite
  # Expected value = TRUE
  
  improper_sol[8] <- any(eigen(fit$psi)$values<0)




  
  names(improper_sol) <- c("reliability_coef",
                          "rho_jh",
                          "p_implied",
                          "theta_jh",
                          "std_lambda",
                          "sigma_implied",
                          "r2",
                          "psi")

  # Improper solution 9 : is p_tilde positive definite
  # Expected value = TRUE

  if (!is.null(fit$p_tilde)){
    improper_sol[9] <- any(eigen(fit$p_tilde)$values<0)
    names(improper_sol)[9] <- "p_tilde"
  }


  
  return(improper_sol= improper_sol)
}