calc_sd_ML <- function(Sigma, x, S_bias, J, R_try, vec_indicators_per_bloc, li_gr, li_non_null_beta_gamma, N, li_C, is_endogenes, grad_eps = 1e-6) {
    inv_Sigma <- solve(Sigma)
    li_deriv_Sigma <- list()
    li_deriv_H <- list()
    for (i in 1:length(x)) {
        x_plus <- x
        x_moins <- x
        x_plus[i] <- x_plus[i] + grad_eps
        x_moins[i] <- x_moins[i] - grad_eps

        li_plus <- reconstruct_blocs(x_plus, S_bias, J, R_try, vec_indicators_per_bloc, li_gr, li_non_null_beta_gamma)
        li_moins <- reconstruct_blocs(x_moins, S_bias, J, R_try, vec_indicators_per_bloc, li_gr, li_non_null_beta_gamma)


        li_Lambda_plus <- li_plus$li_Lambda
        li_phi_exo_plus <- li_plus$li_phi_exo
        li_phi_endo_plus <- li_plus$li_phi_endo
        vec_theta_plus <- li_plus$vec_theta
        li_gamma_plus <- li_plus$li_gamma
        li_beta_plus <- li_plus$li_beta
        is_dag_plus <- li_plus$is_dag
        li_correl_same_diag_plus <- li_plus$li_correl_same_diag


        li_Lambda_moins <- li_moins$li_Lambda
        li_phi_exo_moins <- li_moins$li_phi_exo
        li_phi_endo_moins <- li_moins$li_phi_endo
        vec_theta_moins <- li_moins$vec_theta
        li_gamma_moins <- li_moins$li_gamma
        li_beta_moins <- li_moins$li_beta
        is_dag_moins <- li_moins$is_dag
        li_correl_same_diag_moins <- li_moins$li_correl_same_diag


        Sigma_plus <- calculate_Sigma(li_Lambda_plus, li_phi_exo_plus, li_phi_endo_plus, vec_theta_plus, li_gamma_plus, li_beta_plus, is_dag_plus, li_correl_same_diag_plus, li_C, is_endogenes)
        Sigma_moins <- calculate_Sigma(li_Lambda_moins, li_phi_exo_moins, li_phi_endo_moins, vec_theta_moins, li_gamma_moins, li_beta_moins, is_dag_moins, li_correl_same_diag_moins, li_C, is_endogenes)


        deriv_Sigma <- (Sigma_plus - Sigma_moins) / (2 * grad_eps)
        li_deriv_Sigma[[i]] <- deriv_Sigma


        if (ortho_solve_ML) {
            heq_plus <- heq1(x_plus, S_bias, J, R_try, vec_indicators_per_bloc, li_gr, li_non_null_beta_gamma)
            heq_moins <- heq1(x_moins, S_bias, J, R_try, vec_indicators_per_bloc, li_gr, li_non_null_beta_gamma)
            deriv_H <- (heq_plus - heq_moins) / (2 * grad_eps)
            li_deriv_H[[i]] <- deriv_H
        }
    }
    Info <- matrix(NA, nrow = length(x), ncol = length(x))
    for (i in 1:length(x)) {
        for (j in 1:length(x)) {
            Info[i, j] <- 0.5 * sum(diag(as.matrix(inv_Sigma %*% li_deriv_Sigma[[i]] %*% inv_Sigma %*% li_deriv_Sigma[[j]])))
        }
    }
    Info_avant <- Info

    if (ortho_solve_ML) {
        deriv_H_mat <- do.call(rbind, li_deriv_H)
        I_11 <- Info + deriv_H_mat %*% t(deriv_H_mat)
        I_12 <- deriv_H_mat
        I_21 <- t(deriv_H_mat)
        I_22 <- matrix(0, nrow = nrow(I_21), ncol = ncol(I_12))
        I <- rbind(cbind(I_11, I_12), cbind(I_21, I_22))
        Info <- solve(I)[1:length(x), 1:length(x)]
    }
    cov_mat_uncertainty <- Info / N
    vec_uncertainty <- sqrt(diag(cov_mat_uncertainty))
    reconstructed_uncertainty <- reconstruct_blocs(x, S_bias, J, R_try, vec_indicators_per_bloc, li_gr, li_non_null_beta_gamma)

    return(list(cov_mat_uncertainty = cov_mat_uncertainty, vec_uncertainty = vec_uncertainty, reconstructed_uncertainty = reconstructed_uncertainty))
}
