
#' Detect Improper Solutions in SEM Estimation
#'
#' Checks for nine types of inadmissible or improper solutions in structural
#' equation models, including out-of-range estimates and non-positive definite
#' matrices.
#'
#' @param fit List containing model fit results from `lvm_*()` functions, including:
#'   \itemize{
#'     \item `reliability_coef`: Reliability coefficients for each block
#'     \item `P_IMPLIED`: Implied correlation matrix of latent variables
#'     \item `Ptilde`: Corrected correlation matrix of latent variables
#'     \item `residual_variance`: List of residual variances
#'     \item `std_lambda`: Standardized loadings
#'     \item `SIGMA_IMPLIED`: Implied covariance matrix of observed variables
#'     \item `R2`: R-squared values for endogenous latent variables
#'     \item `psi`: Residual covariance matrix of latent variables
#'   }
#'
#' @return Named logical vector of length 9 indicating presence of each type
#'   of improper solution:
#'   \item{RELIABILITY_COEF}{`TRUE` if any reliability coefficient is outside (0, 1).}
#'   \item{RHO_JH}{`TRUE` if any correlation in P_IMPLIED is outside (-1, 1).}
#'   \item{P_TILDE}{`TRUE` if P_TILDE has negative eigenvalues (not positive definite).}
#'   \item{P_IMPLIED}{`TRUE` if P_IMPLIED has negative eigenvalues (not positive definite).}
#'   \item{THETA_JH}{`TRUE` if any residual variance is negative.}
#'   \item{STD_LAMBDA}{`TRUE` if any standardized loading is outside (-1, 1).}
#'   \item{SIGMA_IMPLIED}{`TRUE` if SIGMA_IMPLIED has negative eigenvalues (not positive definite).}
#'   \item{R2}{`TRUE` if any R-squared is outside (0, 1).}
#'   \item{PSI}{`TRUE` if PSI has negative eigenvalues (not positive definite).}
#'
#' @details
#' A value of `TRUE` indicates an improper solution for that criterion. An
#' admissible solution should have all values set to `FALSE`. Common causes
#' of improper solutions include poor model specification, insufficient sample
#' size, or convergence to a boundary solution.
#'
#' @export
improper <- function(semfc){

  fit <- semfc$estimate
  
  improper_sol <- rep(FALSE, 8)
  
  # Improper solution 1 : number of out of range reliability coefficient
  # Expected value : 0 <=rho_jk <=1

  improper_sol[1] <- any(fit$reliability_coef < 0 |
                        fit$reliability_coef > 1)
  
  
  # Improper solution 2 : number of out of range rho_jk 
  # Expected value : -1 <=rho_jk <=1
  
  improper_sol[2] <- any(abs(fit$P_IMPLIED[upper.tri(fit$P_IMPLIED)])>1)
  

  
  # Improper solution 3 : is P_IMPLIED positive definite
  # Expected value = TRUE
  
  improper_sol[3] <- any(eigen(fit$P_IMPLIED)$values<0)
  
  # Improper solution 4 : is variance of the residual positive
  # Expected value = TRUE
  
  improper_sol[4] <- any(unlist(fit$residual_variance)<0)
  
  # Improper solution 5 : number of out of range cor(y_jh, eta_j)
  # Expected value : -1 <= cor(y_jh, eta_j) <=1
  
  improper_sol[5] <- any(abs(unlist(fit$std_lambda))>1)
  
  # Improper solution 6 : is SIGMA_IMPLIED positive definite
  # Expected value = TRUE
  
  improper_sol[6] <- any(eigen(fit$SIGMA_IMPLIED)$values<0)
  
  # Improper solution 7 : number of out of range R^2_i
  # Expected value : 0 <= R^2_i <=1
  
  improper_sol[7] <- any(fit$R2< 0 | fit$R2>1)
  
  # Improper solution 8 : is PSI positive definite
  # Expected value = TRUE
  
  improper_sol[8] <- any(eigen(fit$psi)$values<0)




  
  names(improper_sol) <- c("RELIABILITY_COEF",
                          "RHO_JH",
                          "P_IMPLIED", 
                          "THETA_JH", 
                          "STD_LAMBDA", 
                          "SIGMA_IMPLIED",
                          "R2",
                          "PSI")

  # Improper solution 9 : is P_TILDE positive definite
  # Expected value = TRUE

  if (!is.null(fit$Ptilde)){
    improper_sol[9] <- any(eigen(fit$Ptilde)$values<0)
    names(improper_sol)[9] <- "P_TILDE"
  }


  
  return(improper_sol= improper_sol)
}