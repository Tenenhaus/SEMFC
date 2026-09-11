calc_P_ij_tilde <- function(S_ij, li_Lambda_star, li_Lambda, li_vec_norm, i, j, R_try) {
    svd_fit <- svd(S_ij, nu = R_try, nv = R_try)
    U <- svd_fit$u
    V <- svd_fit$v
    d <- svd_fit$d[1:R_try]
    lambda_i_star <- li_Lambda_star[[i]]
    lambda_j_star <- li_Lambda_star[[j]]
    mat_ps_gauche <- t(li_Lambda_star[[i]]) %*% U
    mat_ps_droite <- t(li_Lambda_star[[j]]) %*% V
    mat_ps <- mat_ps_gauche + mat_ps_droite
    vec_reorder <- as.integer(clue::solve_LSAP(abs(mat_ps), maximum = TRUE))
    mat_ps_reorder <- mat_ps[, vec_reorder]
    d_good_order <- d[vec_reorder]
    vec_norm_i <- li_vec_norm[[i]]
    vec_norm_j <- li_vec_norm[[j]]
    P_ij_tilde <- diag(d_good_order / (vec_norm_i * vec_norm_j), nrow = R_try, ncol = R_try)
    return(P_ij_tilde)
}
