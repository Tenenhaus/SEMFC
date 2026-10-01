calculate_Sigma <- function(J, R_try, li_Lambda, li_phi_exo, li_phi_endo, vec_theta, li_gamma, li_beta, is_dag, li_correl_same_diag, li_C) {
    # print("go")
    # print(li_Lambda)

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
    return(list(Sigma_implied = Sigma_implied, P_induced = P_induced))
}
