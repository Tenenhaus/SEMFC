gof_ML <- function(x, S, J, R_try, vec_indicators_per_bloc, li_gr, li_non_null_beta_gamma, N, n_endogenes, is_dag) {
    stat <- N * F1(x, S, J, R_try, vec_indicators_per_bloc, li_gr, li_non_null_beta_gamma)
    p <- ncol(S)
    df <- round(p * (p + 1) / 2) - length(x)
    if (ortho_solve_ML) {
        df <- df + round((R_try - 1) * R_try / 2) * J
    }
    p_value <- 1 - stats::pchisq(stat, df)
    return(p_value)
}
