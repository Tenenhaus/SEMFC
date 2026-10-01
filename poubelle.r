local_steps <- function(step_fn, x0, r = 0.02, tol = 1e-7, max_outer = 1000) {
    x <- x0
    for (k in seq_len(max_outer)) {
        fited <- step_fn(x, x - r, x + r)
        x_new <- fited$sol
        print(paste("Iteration", k, "value:", fited$value, "convergence code:", fited$convergence))
        if (max(abs(x_new - x)) < tol) {
            print("optimal attained")
            break # optimum local atteint
        }
        x <- x_new
    }

    if (k == max_outer) {
        print("Maximum number of outer iterations reached without convergence.")
    }
    return(list(sol = x, value = fited$value, convergence = fited$convergence))
}

step_solnp <- function(x, lo, up) {
    n_zeros_eq <- J * R_try * (R_try - 1) / 2
    n_zeros_ineq <- sum(vec_indicators_per_bloc)
    fited <- solnp(
        pars = x, fun = F1, eqfun = heq1, eqB = rep(0, n_zeros_eq),
        LB = pmax(lo, lb), UB = up,
        S = S, J = J, R_try = R_try, vec_indicators_per_bloc = vec_indicators_per_bloc,
        li_gr = fit.svd$li_gr, li_non_null_beta_gamma = li_non_null_beta_gamma, li_C = li_C,
        control = list(trace = 0, tol = 1e-8, delta = 1e-7, rho = 1)
    )
    sol <- fited$pars
    value <- fited$values[length(fited$values)]
    convergence <- fited$convergence
    return(list(sol = sol, value = value, convergence = convergence))
}

index_theta <- (length(big_vec_Lambda) + length(big_vec_phi_exo) + length(big_vec_phi_endo) + 1):(length(big_vec_Lambda) + length(big_vec_phi_exo) + length(big_vec_phi_endo) + length(theta))

lb <- rep(-Inf, length(init_ml_with_S))
lb[index_theta] <- 0

heps <- 1e-7
F1c <- function(p) {
    F1(p,
        S = S, J = J, R_try = R_try,
        vec_indicators_per_bloc = vec_indicators_per_bloc,
        li_gr = fit.svd$li_gr,
        li_non_null_beta_gamma = li_non_null_beta_gamma,
        li_C = li_C
    )
}

heq1c <- function(p) {
    heq1(p,
        S = S, J = J, R_try = R_try,
        vec_indicators_per_bloc = vec_indicators_per_bloc,
        li_gr = fit.svd$li_gr,
        li_non_null_beta_gamma = li_non_null_beta_gamma,
        li_C = li_C
    )
}

grad_F1c <- function(p) nloptr::nl.grad(p, F1c, heps = heps)
jac_heq1c <- function(p) nloptr::nl.jacobian(p, heq1c, heps = heps)

step_nloptr <- function(x, lo, up) {
    fited <- nloptr::nloptr(
        x0 = x, eval_f = F1c, eval_grad_f = grad_F1c,
        lb = pmax(lo, lb), ub = up,
        eval_g_eq = heq1c, eval_jac_g_eq = jac_heq1c,
        opts = list(algorithm = "NLOPT_LD_SLSQP", xtol_rel = 1e-8)
    )
    sol <- fited$solution
    value <- fited$objective
    convergence <- fited$status
    return(list(sol = sol, value = value, convergence = convergence))
}

# fit.ml <- local_steps(step_solnp, init_ml_with_S, r = 0.01)
fit.ml <- local_steps(step_nloptr, init_ml_with_S, r = 0.02)

# OTHER METHOD

fit.ml <- nloptr::nloptr(
    x0 = init_ml_with_S,
    eval_f = F1c,
    eval_grad_f = grad_F1c,
    lb = lb,
    eval_g_eq = heq1c,
    eval_jac_g_eq = jac_heq1c,
    opts = list(
        algorithm = "NLOPT_LD_SLSQP", xtol_rel = 1e-8,
        maxeval = 1000, print_level = 1
    )
)

opt_param <- fit.ml$solution
opt_vrais <- fit.ml$objective
convergence <- fit.ml$status
print(paste("convergence code:", convergence))
