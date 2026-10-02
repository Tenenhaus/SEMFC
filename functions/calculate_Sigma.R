calculate_Sigma <- function(J, R_try, li_Lambda, li_phi_exo, li_phi_endo, vec_theta, li_gamma, li_beta, is_dag, li_C, P_implied) {
    li_which_exo_endo <- lapply(li_C, function(C) {
        out <- ind_exo_endo(C)
        return(out)
    })
    Lambda_diag <- Matrix::bdiag(li_Lambda)
    Sigma_implied <- Lambda_diag %*% P_implied %*% t(Lambda_diag)
    diag(Sigma_implied) <- diag(Sigma_implied) + vec_theta
    return(list(Sigma_implied = Sigma_implied, P_implied = P_implied))
}
