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
  li_correl_same_diag <- reconstructed$li_correl_same_diag

  # print("go")
  # print(li_phi_exo)

  calculated <- calculate_Sigma(J, R_try, li_Lambda, li_phi_exo, li_phi_endo, vec_theta, li_gamma, li_beta, is_dag, li_correl_same_diag, li_C)
  Sigma_implied <- calculated$Sigma_implied
  P_induced <- calculated$P_induced
  eigenvalues <- eigen(Sigma_implied, symmetric = TRUE)$values
  # print(paste("min eigenvalue Sigma:", min(eigenvalues)))

  # print(eigen(Sigma_implied)$values, symmetric = TRUE)

  log_vrais <- log(abs(det(Sigma_implied))) + sum(diag(solve(Sigma_implied) %*% S)) -
    log(det(S)) -
    NCOL(S) # - 1e10 * min(min(eigen(Sigma_implied)$values, symmetric = TRUE), 0)


  # print(paste("log_vrais", log_vrais))
  # print(log(det(Sigma_implied)))
  # print(sum(diag(solve(Sigma_implied) %*% S)))

  # print(eigen(Sigma_implied, symmetric = TRUE)$values)

  # if (is.nan(log_vrais)) {
  #   print("ouille")
  #   print(eigen(Sigma_implied, symmetric = TRUE)$values)
  #   stop()
  # }

  # print(li_Lambda[[1]])
  # print(paste("log_vrais", log_vrais))
  if (log_vrais < 0 | min(eigenvalues) < 0) {
    # print("log_vrais < 0")
    stop("erreur")
  }

  # print(round(Lambda_TRUE %*% as.matrix(P_TRUE), 3))

  # print(sum(diag(solve(Sigma_implied) %*% S)))
  # print(ncol(S))

  # print(round(Sigma_implied - S, 3))

  # print(paste("log_vrais", log_vrais))


  #   l1 <- x[1:3]
  #   l2 <- x[4:6]
  #   l3 <- x[7:9]
  #   l4 <- x[10:12]
  #   l5 <- x[13:15]
  #   l6 <- x[16:18]
  #   L <- bdiag(list(l1, l2, l3, l4, l5, l6))
  #   # correlations between exogeneous
  #   P_EXO <- matrix(c(
  #     1, x[19], x[20], x[22],
  #     x[19], 1, x[21], x[23],
  #     x[20], x[21], 1, x[24],
  #     x[22], x[23], x[24], 1
  #   ), 4, 4)
  #   # path coefficients (Gamma en fait)
  #   G <- matrix(c(
  #     x[25], x[26], 0, 0,
  #     0, 0, x[27], x[28]
  #   ), 2, 4, byrow = TRUE)

  #   # path coefficients
  #   B <- matrix(c(
  #     0, x[29],
  #     x[30], 0
  #   ), 2, 2, byrow = TRUE)

  #   # Correlations between endogeneous
  #   P_ENDO <- matrix(c(
  #     1, x[31],
  #     x[31], 1
  #   ), 2, 2)

  #   # Correlations between Latent/emergent variables
  #   R <- rbind(
  #     cbind(P_EXO, P_EXO %*% t(G) %*% t(solve(diag(2) - B))),
  #     cbind(solve(diag(2) - B) %*% G %*% P_EXO, P_ENDO)
  #   )

  #   # cov between MVs
  #   S1 <- matrix(c(
  #     x[32], x[33], x[35],
  #     x[33], x[34], x[36],
  #     x[35], x[36], x[37]
  #   ), 3, 3)
  #   # cov between MVs
  #   S2 <- matrix(c(
  #     x[38], x[39], x[41],
  #     x[39], x[40], x[42],
  #     x[41], x[42], x[43]
  #   ), 3, 3)
  #   # cov between MVs
  #   S3 <- matrix(c(
  #     x[44], x[45], x[47],
  #     x[45], x[46], x[48],
  #     x[47], x[48], x[49]
  #   ), 3, 3)
  #   # cov between MVs
  #   S4 <- matrix(c(
  #     x[50], x[51], x[53],
  #     x[51], x[52], x[54],
  #     x[53], x[54], x[55]
  #   ), 3, 3)

  #   T5 <- diag(x[56:58])
  #   T6 <- diag(x[59:61])

  #   implied_S <- L %*% R %*% t(L) +
  #     bdiag(list(
  #       S1 - l1 %*% t(l1),
  #       S2 - l2 %*% t(l2),
  #       S3 - l3 %*% t(l3),
  #       S4 - l4 %*% t(l4),
  #       T5,
  #       T6
  #     ))


  return(log_vrais)
}
