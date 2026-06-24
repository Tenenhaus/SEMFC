




jac_vech_Sigma <- function(J_Sigma_Lambda, J_Lambda_theta,
                           J_Sigma_R, J_R_theta,
                           J_Theta_theta){
  J <- J_Sigma_Lambda %*% J_Lambda_theta +
    J_Sigma_R %*% J_R_theta + J_Theta_theta

  return(J)
}


jac_Sigma_Lambda <- function(Lambda, R, D_p){

  p <- nrow(Lambda)
  I <- Diagonal(p)

  D_plus <- solve(crossprod(D_p), t(D_p))
  J <- 2 * (D_plus %*% kronecker(Lambda %*% R, I))

  return(J)
}

jac_Sigma_R <- function(Lambda, L_p, D_bar_J){

  # 4. Calcul de la Jacobienne (Équation 7)
  J_Sigma_R <- L_p %*% kronecker(Lambda, Lambda) %*% D_bar_J

  return(J_Sigma_R)
}








####################################################  Jac Lambda / theta #############################################




calcul_J_Lambda_theta <- function(length_other_param, list_Pj) {

  deriv_lambda <- bdiag(list_Pj)
  n_row <- nrow(deriv_lambda)
  block_zeros <- sparseMatrix(i = integer(0),
                             j = integer(0),
                             dims = c(n_row, length_other_param))

  #  [ deriv_lambda | 0 ]
  J_Lambda_theta <- cbind(deriv_lambda, block_zeros)

  return(J_Lambda_theta)
}




####################################################  Jac Theta / theta #############################################

calcul_Jacobienne_Theta_j_reflective <- function(j, lengths_theta, block_sizes, lengths_values_cov) {
  # lengths_theta = c(18,6,4,2,1,30)
  # 1. Calcul de l'opérateur local
  p_j <- block_sizes[j]
  D_local <- calcul_D_vech_diag(p_j)
  n_row <- nrow(D_local)


  n_theta_other <- tail(cumsum(lengths_theta), 2)[1]
  index_debut <- n_theta_other + cumsum(c(0,lengths_values_cov))[j] + 1
  index_fin <- index_debut + p_j - 1

  # 3. Création de la matrice globale vide (que des zéros, format creux)
  J_Globale <- sparseMatrix(i = integer(0),
                            j = integer(0),
                            x = numeric(0),
                            dims = c(n_row, sum(lengths_theta)))

  # 4. Insertion du bloc D à la position exacte
  J_Globale[, index_debut:index_fin] <- D_local

  return(J_Globale)
}

calcul_Jacobienne_Theta_j_formatif <- function(j, lambda_j, block_sizes, lengths_theta, lengths_values_cov) {

  p_j <- block_sizes[j]
  n_row <- lengths_values_cov[j]


  # Bloc Lambda : -2D^+ (lambda_j ⊗ I_pj)
  D <- duplication_matrix(p_j)
  D_plus <- solve(crossprod(D), t(D))
  I_pj <- Diagonal(p_j)

  mat_lambda_sparse <- -2 * D_plus %*% kronecker(lambda_j, I_pj)


  # Bloc Theta : Matrice identité I_{vech_pj}
  mat_theta_sparse <- Diagonal(n_row)


  idx_lam_start <- cumsum(c(0,block_sizes))[j] + 1
  idx_lam_end   <- idx_lam_start + p_j - 1



  offset_theta_global <- tail(cumsum(lengths_theta), 2)[1]
  idx_theta_start <-  offset_theta_global + cumsum(c(0,lengths_values_cov))[j] + 1
  idx_theta_end   <- idx_theta_start + lengths_values_cov[j] - 1


  J_Globale <- sparseMatrix(i = integer(0),
                            j = integer(0),
                            x = numeric(0),
                            dims = c(n_row, sum(lengths_theta)))

  # Insertion de nos deux blocs aux emplacements exacts
  J_Globale[, idx_lam_start:idx_lam_end] <- mat_lambda_sparse
  J_Globale[, idx_theta_start:idx_theta_end] <- mat_theta_sparse

  return(J_Globale)

}


jac_Theta_j_theta <- function(j, lambda_j, block_sizes, lengths_theta, mode, P_j, D_pj) {


  lengths_values_cov <- block_sizes
  lengths_values_cov[mode == "formative"] <- (block_sizes[mode == "formative"]^2 + block_sizes[mode == "formative"]) / 2

  #  (P_j ⊗ P_j) %*% D_pj
  mat_j <- kronecker(P_j, P_j) %*% D_pj

  if (mode[j] == "formative") {
    return(mat_j %*% calcul_Jacobienne_Theta_j_formatif(j,lambda_j, block_sizes, lengths_theta, lengths_values_cov))
  } else if (mode[j] == "reflective") {
    return(mat_j %*% calcul_Jacobienne_Theta_j_reflective(j, lengths_theta, block_sizes, lengths_values_cov))
  }
}



jac_Theta_theta <- function(lambda, block_sizes, lengths_theta, mode, L_p, list_Pj) {
  return(
    L_p %*% Reduce("+",
                   Map(function(j) {
                     jac_Theta_j_theta(
                       j, as.vector(lambda[[j]]), block_sizes, lengths_theta, mode,
                       list_Pj[[j]], duplication_matrix(block_sizes[j])
                     )
                   }, seq_along(block_sizes))
    )
  )
}




####################################################  Jac R / theta #############################################



######################################## Modèle Structurel Non-Récursif #########################





calcul_J_rhoR_endo <- function(P_endo, M_endo, L_bar) {

  # Produit de Kronecker
  Kron_P <- kronecker(P_endo, P_endo)

  # Produit en chaîne
  J_endo <- L_bar %*% Kron_P %*% M_endo

  return(J_endo)
}

calcul_J_rhoR_exo <- function(P_endo, P_exo, C, Gamma, M_exo, L_bar, K_m) {




  I_m2 <- Diagonal(nrow(K_m))
  J_exo <- L_bar %*% (kronecker(P_exo, P_exo) + (I_m2 + K_m) %*% kronecker(P_endo %*% C %*% Gamma, P_exo)) %*% M_exo

  return(J_exo)
}






calcul_J_rhoR_Gamma <- function(P_endo, P_exo, C, Phi_exo, M_Gamma) {

  D_bar <- correlation_duplication_matrix(nrow(P_endo))
  D_bar_plus <- solve(crossprod(D_bar), t(D_bar))



  Kron_prod <- kronecker(P_exo %*% Phi_exo,  P_endo %*% C)

  # Produit en chaîne final
  J_Gamma <- 2 * (D_bar_plus %*% Kron_prod %*% M_Gamma)

  return(J_Gamma)
}

calcul_J_rhoR_B <- function(P_endo, P_exo, C, Gamma, Phi_exo, M_B) {

  D_bar <- correlation_duplication_matrix(nrow(P_endo))
  D_bar_plus <- solve(crossprod(D_bar), t(D_bar))



  Kron_prod <- kronecker(P_exo %*% Phi_exo %*% t(Gamma) %*% t(C), P_endo %*% C)

  # Produit en chaîne final
  J_B <- 2 * (D_bar_plus %*% Kron_prod %*% M_B)

  return(J_B)
}




calcul_J_rhoR_theta <- function(J_exo, J_Gamma, J_B, J_endo, n_lambda, n_Theta) {


  n_row <- nrow(J_endo)


  zeros_lambda <- sparseMatrix(i = integer(0),
                               j = integer(0),
                               x = numeric(0),
                               dims = c(n_row, n_lambda))

  zeros_Theta  <- sparseMatrix(i = integer(0),
                               j = integer(0),
                               x = numeric(0),
                               dims = c(n_row, n_Theta))

  J_rhoR_theta <- cbind(zeros_lambda, J_exo, J_Gamma, J_B, J_endo, zeros_Theta)

  return(J_rhoR_theta)
}









######################################## F #########################



calcul_jac_F_Sigma <- function(Sigma, S) {


  p <- nrow(Sigma)

  Sigma_inv <- solve(Sigma)

  vec_W <- as.vector(Sigma_inv %*% (Sigma - S) %*% Sigma_inv)



  gradient <- t(vec_W) %*% duplication_matrix(p)

  return(gradient)
}


calcul_Gradient <- function(J_F_Sigma, J_Sigma_theta) {



  return(J_F_Sigma %*% J_Sigma_theta)
}






calcul_Hessian <- function(Sigma, D, J_vech) {


  p <- nrow(Sigma)
  W <- solve(Sigma)
  J_vec <- D %*% J_vech
  I_p <- Diagonal(p)

  M1 <-  kronecker(W, I_p) %*% J_vec
  M2 <- kronecker(I_p, W) %*% J_vec


  H <- crossprod(M1, M2)

  return(H)
}




