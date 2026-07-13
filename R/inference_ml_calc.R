




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

calcul_Jacobienne_Theta_j_reflective <- function(p_j, n_total_cols, start_theta_j) {

  # Création du bloc local
  D_local <- calcul_D_vech_diag(p_j)
  n_row <- nrow(D_local)

  index_fin <- start_theta_j + p_j - 1

  # Création globale
  J_Globale <- sparseMatrix(i = integer(0),
                            j = integer(0),
                            x = numeric(0),
                            dims = c(n_row, n_total_cols))
  J_Globale[, start_theta_j:index_fin] <- D_local

  return(J_Globale)
}

calcul_Jacobienne_Theta_j_formatif <- function(lambda_j, p_j, n_row, D_pj, n_total_cols, start_lam_j, start_theta_j) {


  DtD <- crossprod(D_pj)
  D_plus <- solve(DtD, t(D_pj))

  I_pj <- Diagonal(p_j)

  # 2. Blocs locaux
  J_lambda <- -2 * (D_plus %*% kronecker(lambda_j, I_pj))
  J_cov <- Diagonal(n_row)

  # 3. Indices de fin
  idx_lam_end   <- start_lam_j + p_j - 1
  idx_theta_end <- start_theta_j + n_row - 1

  # 4. Création et insertion
  J_Globale <- sparseMatrix(i = integer(0), j = integer(0), x = numeric(0), dims = c(n_row, n_total_cols))
  J_Globale[, start_lam_j:idx_lam_end] <- J_lambda
  J_Globale[, start_theta_j:idx_theta_end] <- J_cov

  return(J_Globale)
}



jac_Theta_j_theta <- function(lambda_j, p_j, n_row, mode_j, P_j, D_pj, n_total_cols, start_lam_j, start_theta_j) {


  mat_j <- kronecker(P_j, P_j) %*% D_pj

  # Routage direct avec les constantes
  if (mode_j == "formative") {
    J_loc <- calcul_Jacobienne_Theta_j_formatif(lambda_j, p_j, n_row, D_pj, n_total_cols, start_lam_j, start_theta_j)
  } else if (mode_j == "reflective") {
    J_loc <- calcul_Jacobienne_Theta_j_reflective(p_j, n_total_cols, start_theta_j)
  }

  return(mat_j %*% J_loc)
}


jac_Theta_theta <- function(lambda, block_sizes, lengths_theta, mode, L_p, list_Pj) {



  n_total_cols <- sum(lengths_theta)
  offset_theta_global <- tail(cumsum(lengths_theta), 2)[1]

  # sizes of the covariance parameters for each block of indicators
  lengths_values_cov <- block_sizes
  idx_form <- mode == "formative"
  lengths_values_cov[idx_form] <- (block_sizes[idx_form]^2 + block_sizes[idx_form]) / 2

  # Pré-calcul EXHAUSTIF de tous les index de départ pour éviter les cumsum dans la boucle
  # On ajoute un 0 au début du cumsum pour avoir l'index "avant" le bloc courant
  starts_lambda <- cumsum(c(0, block_sizes))
  starts_theta  <- offset_theta_global + cumsum(c(0, lengths_values_cov))


  list_matrix <- Map(function(j) {


    jac_Theta_j_theta(
      lambda_j = as.vector(lambda[[j]]),
      p_j = block_sizes[j],
      n_row = lengths_values_cov[j],
      mode_j = mode[j],
      P_j = list_Pj[[j]],
      D_pj =  duplication_matrix(block_sizes[j]),
      n_total_cols = n_total_cols,
      start_lam_j = starts_lambda[j] + 1,
      start_theta_j = starts_theta[j] + 1
    )
  }, seq_along(block_sizes))


  return(L_p %*% Reduce("+", list_matrix))
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


calcul_J_rhoR_Gamma <- function(P_endo, P_exo, C, Phi_exo, M_Gamma, D_bar_plus) {

  Kron_prod <- kronecker(P_exo %*% Phi_exo,  P_endo %*% C)

  # Produit en chaîne final
  J_Gamma <- 2 * (D_bar_plus %*% Kron_prod %*% M_Gamma)

  return(J_Gamma)
}

calcul_J_rhoR_B <- function(P_endo, P_exo, C, Gamma, Phi_exo, M_B, D_bar_plus) {


  Kron_prod <- kronecker(P_exo %*% Phi_exo %*% t(Gamma) %*% t(C), P_endo %*% C)

  # Produit en chaîne final
  J_B <- 2 * (D_bar_plus %*% Kron_prod %*% M_B)

  return(J_B)
}

calcul_J_rhoR_nonrecursif <- function(P_endo, P_exo, C, Gamma, Phi_exo,
                                      M_exo, M_endo, M_Gamma, M_B,
                                      D_bar_plus, L_bar, K_m) {

  J_rhoR_endo <- calcul_J_rhoR_endo(P_endo, M_endo, L_bar)
  J_rhoR_exo <- calcul_J_rhoR_exo(P_endo, P_exo, C, Gamma, M_exo, L_bar, K_m)
  J_rhoR_gamma <- calcul_J_rhoR_Gamma(P_endo, P_exo, C, Phi_exo, M_Gamma, D_bar_plus)
  J_rhoR_B <- calcul_J_rhoR_B(P_endo, P_exo, C, Gamma, Phi_exo, M_B, D_bar_plus)

  return(list(J_rhoR_endo = J_rhoR_endo,
              J_rhoR_exo = J_rhoR_exo,
              J_rhoR_gamma = J_rhoR_gamma,
              J_rhoR_B = J_rhoR_B))
}


######################################## Modèle Structurel Récursif #########################


calculer_K_psi <- function(L_bar, P_endo, C) {
  PC <- P_endo %*% C
  K_psi <- L_bar %*% KhatriRao(PC, PC)

  return(K_psi)
}

calcul_J_rhoR_exo_recursif <- function(M, C, Gamma, M_Phi_exo, L_bar, K_psi, H_inv) {

  # 1. Terme principal : L_bar(M ⊗ M)
  Term1 <- L_bar %*% kronecker(M, M)

  # 2. Terme de correction récursif : K_psi H^{-1} S_diag (C.Gamma ⊗ C.Gamma)
  C_Gamma <- C %*% Gamma
  # Term2 <- K_psi %*% H_inv %*% S_diag %*% kronecker(C_Gamma, C_Gamma)
  Term2 <- K_psi %*% H_inv %*% t(KhatriRao(t(C_Gamma), t(C_Gamma)))

  # 3. Application de la matrice de sélection M_Phi_exo
  J_exo <- (Term1 - Term2) %*% M_Phi_exo

  return(J_exo)
}


calcul_J_rhoR_Gamma_recursif <- function(M, C, Phi_exo, Gamma, P_endo, M_Gamma, D_bar_plus, K_psi, H_inv) {

  # 1. Terme principal : 2*D_bar_plus (M.Phi_exo ⊗ P_endo.C)
  M_Phi <- M %*% Phi_exo
  PC <- P_endo %*% C
  Term1 <- 2 * D_bar_plus %*% kronecker(M_Phi, PC)

  # 2. Terme de correction récursif : K_psi H^{-1} 2*S_diag (C.Gamma.Phi_exo ⊗ C)
  CGP <- C %*% Gamma %*% Phi_exo
  Term2 <- K_psi %*% H_inv %*% (2 * t(KhatriRao(t(CGP), t(C))))

  # 3. Application de la matrice de sélection M_Gamma
  J_Gamma <- (Term1 - Term2) %*% M_Gamma

  return(J_Gamma)
}


calcul_J_rhoR_B_recursif <- function(R, R_endo, P_endo, C, M_B, D_bar_plus, K_psi, H_inv) {

  # 1. Terme principal : 2*D_bar_plus (R.P_endo ⊗ P_endo.C)
  RP <- R %*% P_endo
  PC <- P_endo %*% C
  Term1 <- 2 * D_bar_plus %*% kronecker(RP, PC)

  # 2. Terme de correction récursif : K_psi H^{-1} 2*S_diag (R_endo ⊗ C)
  Term2 <- K_psi %*% H_inv %*% (2 * t(KhatriRao(t(as(R_endo, "dgCMatrix")), t(C))))

  # 3. Application de la matrice de sélection M_B
  J_B <- (Term1 - Term2) %*% M_B

  return(J_B)
}



calcul_J_rhoR_recursif <- function(C, Gamma, Phi_exo, R, R_endo,
                                   P_endo, P_exo,
                                   M_Phi_exo, M_Gamma, M_B,
                                   D_bar_plus, L_bar) {
  H_inv <- solve(C*C)
  M <- P_exo + P_endo %*% C %*% Gamma

  K_psi <- calculer_K_psi(L_bar, P_endo, C)

  J_rhoR_exo <- calcul_J_rhoR_exo_recursif(M, C, Gamma, M_Phi_exo, L_bar, K_psi, H_inv)
  J_rhoR_gamma <- calcul_J_rhoR_Gamma_recursif(M, C, Phi_exo, Gamma, P_endo, M_Gamma, D_bar_plus, K_psi, H_inv)
  J_rhoR_B <- calcul_J_rhoR_B_recursif(R, R_endo, P_endo, C, M_B, D_bar_plus, K_psi, H_inv)

  return(list(J_rhoR_exo = J_rhoR_exo,
              J_rhoR_gamma = J_rhoR_gamma,
              J_rhoR_B = J_rhoR_B))
}



calcul_J_rhoR <- function(dag, C, Gamma, Phi_exo, R, R_endo,
                          P_endo, P_exo,
                          M_exo, M_endo, M_Gamma, M_B,
                          D_bar_plus, L_bar, K_m) {

  if (dag) {
    J_rhoR <- calcul_J_rhoR_recursif(C, Gamma, Phi_exo, R, R_endo,
                                     P_endo, P_exo,
                                     M_exo, M_Gamma, M_B,
                                     D_bar_plus, L_bar)
  } else {
    J_rhoR <- calcul_J_rhoR_nonrecursif(P_endo, P_exo , C, Gamma, Phi_exo,
                                        M_exo, M_endo, M_Gamma, M_B,
                                        D_bar_plus, L_bar, K_m)
  }

  return(J_rhoR)
}





assemble_J_rhoR_theta <- function(n_lambda, n_Theta,
                                  J_rhoR_exo, J_rhoR_gamma, J_rhoR_B, J_rhoR_endo = NULL) {
  n_row <- nrow(J_rhoR_exo)
  zeros_lambda <- sparseMatrix(i = integer(0),
                               j = integer(0),
                               x = numeric(0),
                               dims = c(n_row, n_lambda))

  zeros_Theta  <- sparseMatrix(i = integer(0),
                               j = integer(0),
                               x = numeric(0),
                               dims = c(n_row, n_Theta))
  J_rhoR_theta <- cbind(zeros_lambda, J_rhoR_exo, J_rhoR_gamma, J_rhoR_B, J_rhoR_endo, zeros_Theta)
  return(J_rhoR_theta)
}



compute_J_rhoR_theta <- function(dag, C, Gamma, Phi_exo, R, R_endo, P_endo, P_exo,
                                 M_exo, M_endo,
                                 M_Gamma, M_B,
                                 D_bar_plus, L_bar, K_m,
                                 n_lambda, n_Theta) {
  list_J_rhoR <- calcul_J_rhoR(dag, C, Gamma, Phi_exo, R, R_endo,
                               P_endo, P_exo,
                               M_exo, M_endo,
                               M_Gamma, M_B,
                               D_bar_plus, L_bar, K_m)

  # On combine les scalaires et la liste des matrices en une seule grande liste d'arguments
  args <- c(list(n_lambda = n_lambda, n_Theta = n_Theta), list_J_rhoR)

  # do.call exécute la fonction en lui injectant toute la liste d'un coup
  J_rhoR_theta <- do.call(assemble_J_rhoR_theta, args)

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


compute_gradient <- function(Sigma, S, J_Sigma_theta) {
  J_F <- calcul_jac_F_Sigma(Sigma, S)
  Gradient_ml <- J_F %*% J_Sigma_theta
  return(Gradient_ml)
}




compute_Hessian <- function(Sigma, J_vech) {


  p <- nrow(Sigma)
  W <- solve(Sigma)
  D <- duplication_matrix(p)
  J_vec <- D %*% J_vech
  I_p <- Diagonal(p)

  M1 <-  kronecker(W, I_p) %*% J_vec
  M2 <- kronecker(I_p, W) %*% J_vec


  H <- crossprod(M1, M2)

  return(H)
}





########################################  Constraint #########################



calcul_dh <- function(lambda_j, Sigma_jj, n_row, p_j, D_pj, n_total_cols, start_lam_j, start_theta_j) {


  # 2. Blocs locaux
  J_lambda <- 2 * t(lambda_j) %*% solve(Sigma_jj)

  v_j <- solve(Sigma_jj)%*%lambda_j
  J_cov <- -t(kronecker(v_j, v_j)) %*% D_pj

  # 3. Indices de fin
  idx_lam_end   <- start_lam_j + p_j - 1
  idx_theta_end <- start_theta_j + n_row - 1

  # 4. Création et insertion
  J_Globale <- sparseMatrix(i = integer(0), j = integer(0), x = numeric(0), dims = c(1, n_total_cols))
  J_Globale[, start_lam_j:idx_lam_end] <- J_lambda
  J_Globale[, start_theta_j:idx_theta_end] <- J_cov

  return(J_Globale)
}


compute_gradient_constraint <- function(lambda, Sigma_cov, block_sizes, lengths_theta, mode) {
  n_total_cols <- sum(lengths_theta)
  offset_theta_global <- tail(cumsum(lengths_theta), 2)[1]

  # sizes of the covariance parameters for each block of indicators
  lengths_values_cov <- block_sizes
  idx_form <- mode == "formative"
  lengths_values_cov[idx_form] <- (block_sizes[idx_form]^2 + block_sizes[idx_form]) / 2

  # Pré-calcul EXHAUSTIF de tous les index de départ pour éviter les cumsum dans la boucle
  # On ajoute un 0 au début du cumsum pour avoir l'index "avant" le bloc courant
  starts_lambda <- cumsum(c(0, block_sizes))
  starts_theta  <- offset_theta_global + cumsum(c(0, lengths_values_cov))


  lambda_form <- lambda[mode == "formative"]
  lengths_values_cov_form <- lengths_values_cov[mode == "formative"]
  block_sizes_form <- block_sizes[mode == "formative"]
  starts_lambda_form <- starts_lambda[mode == "formative"]

  list_matrix <- Map(function(j) {
    calcul_dh(lambda_form[[j]], Sigma_cov[[j]], lengths_values_cov_form[j], block_sizes_form[j],
              duplication_matrix(block_sizes_form[j]), n_total_cols,
              starts_lambda_form[j] + 1, starts_theta[j] + 1)
  }, seq_along(block_sizes_form))

  return(do.call(rbind, list_matrix))


}









