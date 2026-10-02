reconstruct_blocs <- function(x, S, J, R_try, vec_indicators_per_bloc, li_gr, li_non_null_beta_gamma, li_C) {
    cumsum_indicators <- cumsum(c(0, vec_indicators_per_bloc))
    li_Lambda <- lapply(1:J, function(j) {
        index_low <- (cumsum_indicators[j] * R_try) + 1
        index_high <- cumsum_indicators[j + 1] * R_try
        Lambda_j <- matrix(x[index_low:index_high], nrow = vec_indicators_per_bloc[j], ncol = R_try)
        return(Lambda_j)
    })

    end_previous <- cumsum_indicators[J + 1] * R_try

    x_cholesky_P <- x[(end_previous + 1):(end_previous + R_try * J * (R_try * J + 1) / 2)]
    cholesky_P <- matrix(0, nrow = R_try * J, ncol = R_try * J)
    cholesky_P[upper.tri(cholesky_P, diag = TRUE)] <- x_cholesky_P
    P_IMPLIED <- t(cholesky_P) %*% cholesky_P

    end_previous <- end_previous + length(x_cholesky_P)

    n_exogenes <- sum(!is_endogenes)
    li_which_exo_endo <- lapply(li_C, function(C) {
        out <- ind_exo_endo(C)
        return(out)
    })


    li_phi_exo <- lapply(1:R_try, function(r) {
        phi_exo_r <- matrix(0, nrow = n_exogenes, ncol = n_exogenes)
        for (j in 1:n_exogenes) {
            for (i in 1:n_exogenes) {
                i_exo <- li_which_exo_endo[[r]]$ind_exo[i]
                j_exo <- li_which_exo_endo[[r]]$ind_exo[j]
                index_i <- (i_exo - 1) * R_try + r
                index_j <- (j_exo - 1) * R_try + r
                phi_exo_r[i, j] <- P_IMPLIED[index_i, index_j]
            }
        }
        return(phi_exo_r)
    })

    n_endogenes <- sum(is_endogenes)


    li_phi_endo <- lapply(1:R_try, function(r) {
        phi_endo_r <- matrix(0, nrow = n_endogenes, ncol = n_endogenes)
        for (i in 1:n_endogenes) {
            for (j in 1:n_endogenes) {
                i_endo <- li_which_exo_endo[[r]]$ind_endo[i]
                j_endo <- li_which_exo_endo[[r]]$ind_endo[j]
                index_i <- (i_endo - 1) * R_try + r
                index_j <- (j_endo - 1) * R_try + r
                phi_endo_r[i, j] <- P_IMPLIED[index_i, index_j]
            }
        }
        return(phi_endo_r)
    })

    is_dag <- sapply(li_gr, function(gr) igraph::is_dag(gr))

    vec_theta <- x[(end_previous + 1):(end_previous + cumsum_indicators[J + 1])]

    end_previous <- end_previous + length(vec_theta)

    cum_sum_gamma <- cumsum(c(0, sapply(li_non_null_beta_gamma, function(x) sum(x$non_null_gamma))))

    li_gamma <- lapply(1:R_try, function(r) {
        index_low <- 1 + cum_sum_gamma[r] + end_previous
        index_high <- cum_sum_gamma[r + 1] + end_previous
        x_usefull <- x[index_low:index_high]
        gamma <- matrix(0, nrow = n_endogenes, ncol = n_exogenes)
        index_current <- 1
        for (j in 1:ncol(li_non_null_beta_gamma[[r]]$non_null_gamma)) {
            for (i in 1:nrow(li_non_null_beta_gamma[[r]]$non_null_gamma)) {
                if (li_non_null_beta_gamma[[r]]$non_null_gamma[i, j]) {
                    gamma[i, j] <- x_usefull[index_current]
                    index_current <- index_current + 1
                }
            }
        }
        return(gamma)
    })

    end_previous <- end_previous + cum_sum_gamma[R_try + 1]

    cumsum_beta <- cumsum(c(0, sapply(li_non_null_beta_gamma, function(x) sum(x$non_null_beta))))

    li_beta <- lapply(1:R_try, function(r) {
        index_low <- 1 + cumsum_beta[r] + end_previous
        index_high <- cumsum_beta[r + 1] + end_previous
        x_usefull <- x[index_low:index_high]
        beta <- matrix(0, nrow = n_endogenes, ncol = n_endogenes)
        index_current <- 1
        for (j in 1:ncol(li_non_null_beta_gamma[[r]]$non_null_beta)) {
            for (i in 1:nrow(li_non_null_beta_gamma[[r]]$non_null_beta)) {
                if (li_non_null_beta_gamma[[r]]$non_null_beta[i, j]) {
                    beta[i, j] <- x_usefull[index_current]
                    index_current <- index_current + 1
                }
            }
        }
        return(beta)
    })


    return(list(li_Lambda = li_Lambda, li_phi_exo = li_phi_exo, li_phi_endo = li_phi_endo, vec_theta = vec_theta, li_gamma = li_gamma, li_beta = li_beta, is_dag = is_dag, P_implied = P_IMPLIED))
}
