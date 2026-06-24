#' Final Function to Compute J_Sigma_theta
#'
#' This function computes the final Jacobian matrix J_Sigma_theta with minimal inputs.
#' It handles all intermediate calculations including structural and measurement models.
#'
#' @param lambda List of factor loadings matrices for each block
#' @param block_sizes Numeric vector indicating the size of each block
#' @param lengths_theta Numeric vector indicating parameter counts for each theta component
#' @param mode Character vector indicating mode ('reflective' or 'formative') for each block
#' @param R Correlation matrix of latent variables
#' @param rel_matrix Relationship matrix (C matrix) for structural model
#' @param m_exo Number of exogenous variables
#' @param PHI Covariance matrix of exogenous variables
#' @param BETA Regression matrix of structural model
#' @param GAMMA Regression matrix from exogenous to endogenous
#'
#' @return Matrix J_Sigma_theta (Jacobian of vech(Sigma) with respect to theta)
#'
#' @examples
#' \dontrun{
#' J_Sigma_theta_final <- compute_J_Sigma_theta(
#'   lambda = lambda,
#'   block_sizes = c(3,3,3,3,3,3),
#'   lengths_theta = c(18,6,4,2,1,30),
#'   mode = mode,
#'   R = R,
#'   rel_matrix = rel_matrix,
#'   m_exo = 4,
#'   PHI = PHI,
#'   BETA = BETA,
#'   GAMMA = GAMMA
#' )
#' }
#'
#' @export
compute_J_Sigma_theta <- function(block_sizes, lengths_theta, mode, R,
                                   rel_matrix, lambda, PHI, BETA, GAMMA) {


  C <- solve(diag(nrow(BETA)) - BETA)
  which_exo_endo <- ind_exo_endo(rel_matrix)
  ind_exo <- which_exo_endo$ind_exo
  ind_endo <- which_exo_endo$ind_endo
  J <- ncol(rel_matrix)  # Total number of latent variables
  p <- sum(block_sizes)  # Total number of observed variables

  D_p <- duplication_matrix(p)
  D_bar_J <- correlation_duplication_matrix(J)
  D_bar_plus <- solve(crossprod(D_bar_J), t(D_bar_J))
  L_p <- elimination_matrix(p)
  L_bar <- correlation_elimination_matrix(J)



  proj_matrices <- P_exo_endo(ind_exo, ind_endo)
  P_exo <- proj_matrices$P_exo
  P_endo <- proj_matrices$P_endo


  beta_gamma_matrices <- M_beta_gamma(rel_matrix)
  M_gamma <- beta_gamma_matrices$M_gamma
  M_beta <- beta_gamma_matrices$M_beta

  list_Pj <- generate_Pj(block_sizes)
  
  # 1. Lambda Jacobian
  Lambda <- bdiag(lambda)
  J_Sigma_Lambda <- jac_Sigma_Lambda(Lambda, R, D_p)
  J_Lambda_theta <- calcul_J_Lambda_theta(sum(lengths_theta) - p , list_Pj)
  

  
  # 3. Theta Jacobian (measurement model)
  J_Theta_theta <- jac_Theta_theta(lambda, block_sizes, lengths_theta, mode, L_p, list_Pj)
  
  # 2. Correlation Jacobian

  J_Sigma_R <- jac_Sigma_R(Lambda, L_p, D_bar_J)

  J_rhoR_endo <- calcul_J_rhoR_endo(P_endo, 
                                    M_endo = correlation_duplication_matrix(length(ind_endo)),
                                    L_bar = L_bar)
  J_rhoR_exo <- calcul_J_rhoR_exo(P_endo, P_exo, C, GAMMA,
                                  M_exo = correlation_duplication_matrix(length(ind_exo)),
                                  L_bar = L_bar,
                                  K_m = commutation_matrix(J))
  J_rhoR_gamma <- calcul_J_rhoR_Gamma(P_endo, P_exo, C, PHI, M_gamma, D_bar_plus)
  J_rhoR_B <- calcul_J_rhoR_B(P_endo, P_exo, C, GAMMA, PHI, M_beta, D_bar_plus)
  
  # 4e. Combine structural model Jacobians
  J_rhoR_theta <- calcul_J_rhoR_theta(J_rhoR_exo, J_rhoR_gamma, J_rhoR_B, J_rhoR_endo,
                                      n_lambda = p,
                                      n_Theta = tail(lengths_theta, 1))
  
  # 5. Combine all components

  J_Sigma_theta <- jac_vech_Sigma(J_Sigma_Lambda, J_Lambda_theta, J_Sigma_R, J_rhoR_theta, J_Theta_theta)
  
  return(J_Sigma_theta)
}



compute_gradient_hessian <- function(lambda, block_sizes, lengths_theta, mode, R,
                                   rel_matrix, m_exo, PHI, BETA, GAMMA, SIGMA, S) {
  # Compute J_Sigma_theta
  J_Sigma_theta <- compute_J_Sigma_theta(block_sizes, lengths_theta, mode,
                                         R, rel_matrix, m_exo, lambda, PHI, BETA, GAMMA)
  # Compute Gradient
  Gradient <- compute_gradient(Sigma, S, J_Sigma_theta)
  # Compute Hessian
  H <- compute_Hessian(SIGMA, J_Sigma_theta)

  return(list(Gradient = Gradient, Hessian = H))

}
