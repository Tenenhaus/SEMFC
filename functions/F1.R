# ########################################
# # Objective function for ML estimation #
# ########################################

F1 <- function(x, S, J, R_try, vec_indicators_per_bloc, li_gr, li_non_null_beta_gamma, show_eigen = FALSE) {
  #   # loadings

  reconstructed <- reconstruct_blocs(x, S, J, R_try, vec_indicators_per_bloc, li_gr, li_non_null_beta_gamma)

  # li_Lambda = li_Lambda, li_phi_exo = li_phi_exo, li_phi_endo = li_phi_endo, vec_theta = vec_theta, li_gamma = li_gamma, li_beta = li_beta
  li_Lambda <- reconstructed$li_Lambda
  li_phi_exo <- reconstructed$li_phi_exo
  li_phi_endo <- reconstructed$li_phi_endo
  vec_theta <- reconstructed$vec_theta
  li_gamma <- reconstructed$li_gamma
  li_beta <- reconstructed$li_beta
  is_dag <- reconstructed$is_dag
  li_correl_same_diag <- reconstructed$li_correl_same_diag


  Sigma_implied <- calculate_Sigma(li_Lambda, li_phi_exo, li_phi_endo, vec_theta, li_gamma, li_beta, is_dag, li_correl_same_diag, li_C, is_endogenes)

  if (min(eigen(Sigma_implied, symmetric = TRUE)$values) < 1e-8) {
    log_vrais <- Inf
    # print(eigen(P_induced, symmetric = TRUE)$values)
  } else {
    log_vrais <- log(det(Sigma_implied)) + sum(diag(solve(Sigma_implied) %*% S)) -
      log(det(S)) -
      NCOL(S)
  }

  # print(paste("log_vrais", log_vrais))


  if (is.nan(log_vrais)) {
    print("ouille")
    print(eigen(Sigma_implied)$values, symmetric = TRUE)
    stop()
  }

  # print(li_Lambda[[1]])
  # print(paste("log_vrais", log_vrais))
  # print(min(eigen(P_induced)$values, symmetric = TRUE))
  if (show_eigen) {
    diff <- sum(abs(P_induced - t(P_induced)))
    print(diff)
    print(P_induced)
  }

  if (log_vrais < -1e-2) {
    print(min(eigen(Sigma_implied, symmetric = TRUE)$values))
    print(paste("log_vrais", log_vrais))
    print("log_vrais < -1e-2")
    # if (log_vrais < -1) {
    #   stop("erreur")
    # }
  }


  return(log_vrais)
}
