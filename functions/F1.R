# ########################################
# # Objective function for ML estimation #
# ########################################

F1 <- function(x, S, J, R_try, vec_indicators_per_bloc, li_gr, li_non_null_beta_gamma, li_C) {
  #   # loadings

  reconstructed <- reconstruct_blocs(x, S, J, R_try, vec_indicators_per_bloc, li_gr, li_non_null_beta_gamma, li_C)

  # li_Lambda = li_Lambda, li_phi_exo = li_phi_exo, li_phi_endo = li_phi_endo, vec_theta = vec_theta, li_gamma = li_gamma, li_beta = li_beta
  li_Lambda <- reconstructed$li_Lambda
  li_phi_exo <- reconstructed$li_phi_exo
  li_phi_endo <- reconstructed$li_phi_endo
  vec_theta <- reconstructed$vec_theta
  li_gamma <- reconstructed$li_gamma
  li_beta <- reconstructed$li_beta
  is_dag <- reconstructed$is_dag
  P_implied <- reconstructed$P_implied


  # print("go")
  # print(li_phi_exo)

  calculated <- calculate_Sigma(J, R_try, li_Lambda, li_phi_exo, li_phi_endo, vec_theta, li_gamma, li_beta, is_dag, li_C, P_implied)
  Sigma_implied <- calculated$Sigma_implied
  eigenvalues <- eigen(Sigma_implied, symmetric = TRUE)$values
  # print(paste("min eigenvalue Sigma:", min(eigenvalues)))

  # print(eigen(Sigma_implied)$values, symmetric = TRUE)


  log_vrais <- log(det(Sigma_implied)) + sum(diag(solve(Sigma_implied) %*% S)) -
    log(det(S)) - NCOL(S)


  print(paste("log_vrais", log_vrais))
  # if (log_vrais > 10) {
  #   print(Sigma_implied)
  #   stop("log_vrais > 10")
  # }


  # print(log(det(Sigma_implied)))
  # print(sum(diag(solve(Sigma_implied) %*% S)))

  print(paste("min eigenvalue Sigma:", min(eigen(Sigma_implied, symmetric = TRUE)$values)))

  # if (is.nan(log_vrais)) {
  #   print("ouille")
  #   print(eigen(Sigma_implied, symmetric = TRUE)$values)
  #   stop()
  # }

  # print(li_Lambda[[1]])
  # print(paste("log_vrais", log_vrais))
  if (log_vrais < -1e-3) {
    print(eigen(Sigma_implied, symmetric = TRUE)$values)
    stop("erreur")
  }


  return(log_vrais)
}
