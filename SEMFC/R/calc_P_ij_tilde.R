calc_P_ij_tilde <- function(S_ij, li_Lambda_star, li_Lambda, li_vec_norm, i, j, R_try, ortho, S_ii, S_jj) {
    if (ortho) {
        svd_fit <- svd(S_ij, nu = R_try, nv = R_try)
        U <- svd_fit$u
        V <- svd_fit$v
        d <- svd_fit$d[1:R_try]
        lambda_i_star <- li_Lambda_star[[i]]
        lambda_j_star <- li_Lambda_star[[j]]
        mat_ps_gauche <- abs(t(li_Lambda_star[[i]]) %*% U)
        mat_ps_droite <- abs(t(li_Lambda_star[[j]]) %*% V)
        mat_ps <- mat_ps_gauche + mat_ps_droite
        vec_reorder <- as.integer(clue::solve_LSAP(abs(mat_ps), maximum = TRUE))
        mat_ps_reorder <- mat_ps[, vec_reorder]
        # print(mat_ps_reorder)
        d_good_order <- d[vec_reorder]
        vec_norm_i <- li_vec_norm[[i]]
        vec_norm_j <- li_vec_norm[[j]]
        vec_p_ij_tilde <- d_good_order / (vec_norm_i * vec_norm_j)

        D_i_inv <- diag(1 / li_vec_norm[[i]], nrow = R_try, ncol = R_try)
        D_j_inv <- diag(1 / li_vec_norm[[j]], nrow = R_try, ncol = R_try)
        good_signs <- sign(diag(D_i_inv %*% t(li_Lambda_star[[i]]) %*% S_ij %*% li_Lambda_star[[j]] %*% D_j_inv))
        vec_p_ij_tilde <- vec_p_ij_tilde * good_signs

        P_ij_tilde <- diag(vec_p_ij_tilde, nrow = R_try, ncol = R_try)
    } else {
        gram_i <- t(li_Lambda[[i]]) %*% li_Lambda[[i]]
        gram_j <- t(li_Lambda[[j]]) %*% li_Lambda[[j]]
        A <- gram_i * gram_j
        C <- sapply(1:R_try, function(r) {
            lambda_i_r <- li_Lambda[[i]][, r]
            lambda_j_r <- li_Lambda[[j]][, r]
            c_r <- as.vector(t(lambda_i_r) %*% S_ij %*% lambda_j_r)
            return(c_r)
        })
        diag_tilde <- solve(A, C)
        P_ij_tilde <- diag(diag_tilde, nrow = R_try, ncol = R_try)
    }

    return(P_ij_tilde)
}
