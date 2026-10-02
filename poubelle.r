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

### RECALCULATE P USELESS

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
