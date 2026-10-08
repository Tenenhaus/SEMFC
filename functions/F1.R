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

  li_which_exo_endo <- lapply(li_C, function(C) {
    out <- ind_exo_endo(C)
    return(out)
  })


  li_P <- lapply(1:R_try, function(r) {
    H_index <- li_which_exo_endo[[r]]$ind_exo
    J_index <- li_which_exo_endo[[r]]$ind_endo
    phi_exo <- li_phi_exo[[r]]
    gamma <- li_gamma[[r]]
    beta <- li_beta[[r]]
    P_12 <- phi_exo %*% t(gamma) %*% t(solve(diag(NROW(beta)) - beta))
    P_21 <- t(P_12)
    if (is_dag[r]) {
      P_22 <- matrix(0, nrow = NROW(beta), ncol = NROW(beta))
    } else {
      phi_endo <- li_phi_endo[[r]]
      P_22 <- phi_endo
    }
    P_shuffled <- rbind(
      cbind(phi_exo, P_12),
      cbind(P_21, P_22)
    )
    P <- matrix(0, nrow = NCOL(P_shuffled), ncol = NCOL(P_shuffled))
    P[c(H_index, J_index), c(H_index, J_index)] <- P_shuffled
    return(P)
  })

  n_endogenes <- sum(is_endogenes)

  if (any(is_dag)) {
    for (r in 1:R_try) {
      if (is_dag[r]) {
        J_index <- li_which_exo_endo[[r]]$ind_endo
        H_index <- li_which_exo_endo[[r]]$ind_exo
        beta <- li_beta[[r]]
        gamma <- li_gamma[[r]]
        PI <- solve(diag(NROW(beta)) - beta)
        P <- li_P[[r]]
        diag_psi <- sapply(1:n_endogenes, function(i) {
          bg_i <- c(beta[i, ], gamma[i, ])
          psi_ii <- 1 - t(bg_i) %*% P[c(J_index, H_index), c(J_index, H_index)] %*% bg_i
          return(psi_ii)
        })
        mat_psi <- diag(diag_psi, nrow = length(diag_psi), ncol = length(diag_psi))
        P_22 <- PI %*% (gamma %*% P[H_index, H_index] %*% t(gamma) + mat_psi) %*% t(PI)
        P_11 <- P[H_index, H_index]
        P_12 <- P[H_index, J_index]
        P_21 <- P[J_index, H_index]
        P_shuffled <- rbind(
          cbind(P_11, P_12),
          cbind(P_21, P_22)
        )
        P <- matrix(0, nrow = NCOL(P_shuffled), ncol = NCOL(P_shuffled))
        P[c(H_index, J_index), c(H_index, J_index)] <- P_shuffled
        li_P[[r]] <- P
      }
    }
  }

  # OUILLE il manque les corrélations entre P de même rang...?


  P_induced <- matrix(0, nrow = J * R_try, ncol = J * R_try)
  for (i in 1:NCOL(li_P[[1]])) {
    for (j in 1:NCOL(li_P[[1]])) {
      if (i != j) {
        for (r in 1:R_try) {
          P_induced[(i - 1) * R_try + r, (j - 1) * R_try + r] <- li_P[[r]][i, j]
        }
      } else {
        P_induced[((i - 1) * R_try + 1):(i * R_try), ((j - 1) * R_try + 1):(j * R_try)] <- li_correl_same_diag[[i]]
      }
    }
  }


  Lambda_diag <- Matrix::bdiag(li_Lambda)
  Sigma_implied <- Lambda_diag %*% P_induced %*% t(Lambda_diag)
  diag(Sigma_implied) <- diag(Sigma_implied) + vec_theta


  log_vrais <- log(det(Sigma_implied)) + sum(diag(solve(Sigma_implied) %*% S)) -
    log(det(S)) -
    NCOL(S)

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
