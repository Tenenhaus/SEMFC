# restrictions for the minimization of the Loglikelihood function
heq1 <- function(x, S, J, R_try, vec_indicators_per_bloc, li_gr, li_non_null_beta_gamma, li_C) {
  reconstructed <- reconstruct_blocs(x, S, J, R_try, vec_indicators_per_bloc, li_gr, li_non_null_beta_gamma, li_C)


  # Orthogonalité des lambda_i
  li_Lambda <- reconstructed$li_Lambda
  li_orthogonalities <- lapply(1:J, function(j) {
    mat_diag <- t(li_Lambda[[j]]) %*% li_Lambda[[j]]
    orthogonalities <- mat_diag[upper.tri(mat_diag)]
    return(orthogonalities)
  })
  vec_orthogonalities <- do.call("c", li_orthogonalities)

  h <- c(vec_orthogonalities)

  # Elements nuls de P_implied

  P_implied <- reconstructed$P_implied
  is_zero_P <- matrix(TRUE, nrow = R_try * J, ncol = R_try * J)
  for (j in 1:J) {
    for (k in 1:J) {
      if (j != k) {
        for (r in 1:R_try) {
          is_zero_P[(j - 1) * R_try + r, (k - 1) * R_try + r] <- FALSE
        }
      }
      if (j == k) {
        index_low <- (j - 1) * R_try + 1
        index_high <- j * R_try
        is_zero_P[index_low:index_high, index_low:index_high] <- FALSE
      }
    }
  }


  is_zero_P_upper_tri <- is_zero_P[upper.tri(is_zero_P, diag = FALSE)]
  P_upper_tri <- P_implied[upper.tri(P_implied, diag = FALSE)]
  h <- c(h, P_upper_tri[is_zero_P_upper_tri])

  # diagonale - 1 = 0
  vec_diag_null <- diag(P_implied) - 1
  h <- c(h, vec_diag_null)

  ## Lien entre Beta, Gamma et P_implied
  li_beta <- reconstructed$li_beta
  li_gamma <- reconstructed$li_gamma

  li_P_r <- lapply(1:R_try, function(r) {
    P_r <- matrix(0, nrow = J, ncol = J)
    for (i in 1:J) {
      for (j in 1:J) {
        P_r[i, j] <- P_implied[(i - 1) * R_try + r, (j - 1) * R_try + r]
      }
    }
    return(P_r)
  })

  li_which_exo_endo <- lapply(li_C, function(C) {
    out <- ind_exo_endo(C)
    return(out)
  })

  li_vec_constraints_non_diag <- lapply(1:R_try, function(r) {
    beta <- li_beta[[r]]
    gamma <- li_gamma[[r]]
    P_r <- li_P_r[[r]]
    ind_exo <- li_which_exo_endo[[r]]$ind_exo
    ind_endo <- li_which_exo_endo[[r]]$ind_endo
    P_r_non_diag <- P_r[ind_exo, ind_endo]
    delta <- diag(NROW(beta)) - t(beta)
    non_diag_constraints <- P_r_non_diag - P_r[ind_exo, ind_exo] %*% t(gamma) %*% solve(delta)
    return(non_diag_constraints)
  })

  vec_constraints_non_diag <- do.call("c", li_vec_constraints_non_diag)
  h <- c(h, vec_constraints_non_diag)

  ## Contrainte sur les dag

  if (any(sapply(li_gr, function(gr) igraph::is_dag(gr)))) {
    li_vec_constraints_dag <- lapply(1:R_try, function(r) {
      if (igraph::is_dag(li_gr[[r]])) {
        beta <- li_beta[[r]]
        gamma <- li_gamma[[r]]
        P_r <- li_P_r[[r]]
        ind_exo <- li_which_exo_endo[[r]]$ind_exo
        ind_endo <- li_which_exo_endo[[r]]$ind_endo
        delta <- diag(NROW(beta)) - t(beta)
        PSI <- matrix(0, NCOL(beta), NCOL(beta))
        for (i in 1:NCOL(beta)) {
          bg <- c(beta[i, ], gamma[i, ])
          PSI[i, i] <- 1 - as.vector(t(bg) %*% P_r[c(ind_endo, ind_exo), c(ind_endo, ind_exo)] %*% bg)
        }

        constraints_dag <- P_r[ind_endo, ind_endo] - t(delta) %*% (gamma %*% P_r[ind_exo, ind_exo] %*% t(gamma) + PSI) %*% delta
        return(constraints_dag[upper.tri(constraints_dag, diag = FALSE)])
      } else {
        return(c())
      }
    })

    vec_constraints_dag <- do.call("c", li_vec_constraints_dag)
    h <- c(h, vec_constraints_dag)
  }
  print(paste("eqfun max abs:", max(abs(h))))
  ecart <- abs(init_ml_with_S - x)
  # print(paste("eqfun max abs diff:", sum(abs(ecart))))
  # print(eigen(P_implied, symmetric = TRUE)$values)
  if (max(abs(h)) > 1) {
    print("h constraints too large, stopping")
    stop()
  }

  return(h)
}

ineqfun <- function(x, S, J, R_try, vec_indicators_per_bloc, li_gr, li_non_null_beta_gamma, li_C) {
  reconstructed <- reconstruct_blocs(x, S, J, R_try, vec_indicators_per_bloc, li_gr, li_non_null_beta_gamma, li_C)
  li_Lambda <- reconstructed$li_Lambda
  li_phi_exo <- reconstructed$li_phi_exo
  li_phi_endo <- reconstructed$li_phi_endo
  vec_theta <- reconstructed$vec_theta
  li_gamma <- reconstructed$li_gamma
  li_beta <- reconstructed$li_beta
  is_dag <- reconstructed$is_dag
  P_implied <- reconstructed$P_implied

  # calculated <- calculate_Sigma(J, R_try, li_Lambda, li_phi_exo, li_phi_endo, vec_theta, li_gamma, li_beta, is_dag, li_C, P_implied)
  # Sigma_implied <- calculated$Sigma_implied
  # min_values <- min(eigen(Sigma_implied, symmetric = TRUE)$values)
  h <- c(vec_theta)
  return(h)
}
