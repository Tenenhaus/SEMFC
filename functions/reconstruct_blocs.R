reconstruct_blocs <- function(x, S, J, R_try, vec_indicators_per_bloc, li_gr, li_non_null_beta_gamma) {
    cumsum_indicators <- cumsum(c(0, vec_indicators_per_bloc))
    li_Lambda <- lapply(1:J, function(j) {
        index_low <- (cumsum_indicators[j] * R_try) + 1
        index_high <- cumsum_indicators[j + 1] * R_try
        Lambda_j <- matrix(x[index_low:index_high], nrow = vec_indicators_per_bloc[j], ncol = R_try)
        return(Lambda_j)
    })

    end_previous <- cumsum_indicators[J + 1] * R_try

    n_exogenes <- sum(!is_endogenes)

    li_phi_exo <- lapply(1:R_try, function(r) {
        index_low <- 1 + end_previous + (r - 1) * (n_exogenes * (n_exogenes - 1)) %/% 2
        index_high <- index_low + (n_exogenes * (n_exogenes - 1)) %/% 2 - 1
        x_usefull <- x[index_low:index_high]
        phi_exo_r <- matrix(0, nrow = n_exogenes, ncol = n_exogenes)
        phi_exo_r[upper.tri(phi_exo_r)] <- x_usefull
        phi_exo_r <- phi_exo_r + t(phi_exo_r)
        diag(phi_exo_r) <- 1
        return(phi_exo_r)
    })

    n_endogenes <- sum(is_endogenes)

    end_previous <- end_previous + R_try * (n_exogenes * (n_exogenes - 1)) %/% 2

    li_phi_endo <- lapply(1:R_try, function(r) {
        gr <- li_gr[[r]]
        if (!igraph::is_dag(gr)) {
            index_low <- 1 + end_previous + (r - 1) * (n_endogenes * (n_endogenes - 1)) %/% 2
            index_high <- index_low + (n_endogenes * (n_endogenes - 1)) %/% 2 - 1
            x_usefull <- x[index_low:index_high]
            phi_endo_r <- matrix(0, nrow = n_endogenes, ncol = n_endogenes)
            phi_endo_r[upper.tri(phi_endo_r)] <- x_usefull
            phi_endo_r <- phi_endo_r + t(phi_endo_r)
            diag(phi_endo_r) <- 1
        } else {
            phi_endo_r <- NULL
        }
        return(phi_endo_r)
    })

    is_dag <- sapply(li_gr, function(gr) igraph::is_dag(gr))

    end_previous <- end_previous + sum(!is_dag) * (n_endogenes * (n_endogenes - 1)) %/% 2

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

    end_previous <- end_previous + cumsum_beta[R_try + 1]

    li_correl_same_diag <- lapply(1:J, function(j) {
        index_low <- 1 + end_previous + (j - 1) * ((R_try * (R_try - 1)) %/% 2)
        index_high <- index_low + ((R_try * (R_try - 1)) %/% 2) - 1
        x_usefull <- x[index_low:index_high]
        # print(length(x_usefull))
        mat_correl_same_diag <- matrix(0, nrow = R_try, ncol = R_try)
        mat_correl_same_diag[upper.tri(mat_correl_same_diag)] <- x_usefull
        mat_correl_same_diag <- mat_correl_same_diag + t(mat_correl_same_diag)
        diag(mat_correl_same_diag) <- 1
        return(mat_correl_same_diag)
    })


    return(list(li_Lambda = li_Lambda, li_phi_exo = li_phi_exo, li_phi_endo = li_phi_endo, vec_theta = vec_theta, li_gamma = li_gamma, li_beta = li_beta, is_dag = is_dag, li_correl_same_diag = li_correl_same_diag))
}
