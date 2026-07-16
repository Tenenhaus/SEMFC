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
#' @param Phi_exo Covariance matrix of exogenous variables
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
#'   Phi_exo = Phi_exo,
#'   BETA = BETA,
#'   GAMMA = GAMMA
#' )
#' }
#'
#' @export
compute_J_Sigma_theta <- function(block_sizes, lengths_theta, mode, dag, R,
                                   rel_matrix, lambda, BETA, GAMMA) {


  C <- solve(diag(nrow(BETA)) - BETA)
  which_exo_endo <- ind_exo_endo(rel_matrix)
  ind_exo <- which_exo_endo$ind_exo
  ind_endo <- which_exo_endo$ind_endo
  Phi_exo <- R[ind_exo, ind_exo, drop = FALSE]
  Phi_endo <- R[ind_endo, ind_endo, drop = FALSE]
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


  J_rhoR_theta <- compute_J_rhoR_theta(dag, C, GAMMA, Phi_exo, R, Phi_endo,
                                       P_endo, P_exo,
                                       M_exo = correlation_duplication_matrix(length(ind_exo)),
                                       M_endo = correlation_duplication_matrix(length(ind_endo)),
                                       M_gamma, M_beta,
                                       D_bar_plus, L_bar = L_bar, K_m = commutation_matrix(J),
                                       n_lambda = p,
                                       n_Theta = tail(lengths_theta, 1))




  J_Sigma_theta <- jac_vech_Sigma(J_Sigma_Lambda, J_Lambda_theta, J_Sigma_R, J_rhoR_theta, J_Theta_theta)

  return(J_Sigma_theta)
}



grad_F1 <- function (x, S, model) {
  est <- lvm_ml(x, model, jac = FALSE)
  R <- est$p_implied
  lambda <- est$lambda
  BETA <- est$beta
  GAMMA <- est$gamma
  Sigma <- est$sigma_implied

  J_Sigma_theta <- compute_J_Sigma_theta(model$block_sizes, model$lengths_theta, model$mode,model$dag,
                                         R, model$relation_matrix, lambda, BETA, GAMMA)
  # Compute Gradient
  Gradient <- compute_gradient(Sigma, S, J_Sigma_theta)

  return(as.vector(Gradient))
}


hess_F1 <- function(x, S, model) {
  est <- lvm_ml(x, model, jac = FALSE)
  R <- est$p_implied
  lambda <- est$lambda
  BETA <- est$beta
  GAMMA <- est$gamma
  Sigma <- est$sigma_implied

  J_Sigma_theta <- compute_J_Sigma_theta(model$block_sizes, model$lengths_theta, model$mode,model$dag,
                                         R, model$relation_matrix, lambda, BETA, GAMMA)

  H <- compute_Hessian(Sigma, J_Sigma_theta)

  return(as.matrix(H))
}
