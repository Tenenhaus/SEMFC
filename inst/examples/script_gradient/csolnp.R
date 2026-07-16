

devtools::load_all()

library('nloptr')

prob <- solnp_problem_suite(number = 10)
sol <- csolnp(pars = prob$start, fn = prob$fn, gr = prob$gr, eq = prob$eq_fn, eq_b = prob$eq_b,
              eq_jac = prob$eq_jac, ineq_fn = prob$ineq_fn, ineq_lower = prob$ineq_lower,
              ineq_upper = prob$ineq_upper, ineq_jac = prob$ineq_jac, lower = prob$lower,
              upper = prob$upper)
#> Warning in solnp_problem_setup(pars, fn, gr, eq_fn, eq_b, eq_jac, ineq_fn, :
#> lower inequality values violated with initial pars.
print(prob$name)
#> [1] "hs10"
print(c("convergence" = sol$convergence))
#> convergence
#>           0
print(c("csolnp objective" = sol$objective, "best objective" = prob$best_fn))
#> csolnp objective   best objective
#>       -0.9999991       -1.0000000
print(sol$elapsed_time)


source('inst/model/model_mixed.R')
Y <- Y_2
X <- X_2

set.seed(27)
modelsvd <- SemFC$new(data=Y, relation_matrix = C, mode=mode, estimator = "svd")

modelsvd$fit(infer = F)


init <- modelsvd$get_estimate('theta')
model <- modelsvd$get_model()
S <- cov(X)
r = 4
tol = 1e-08

system.time({



result <- solnp(pars = init,
                fun=F1, eqfun=heq1,
                eqB = rep(0,r), S = S, model = model,
                control = list(trace = 0, tol = tol))
})

f <- function (x) {
  return(F1(x, S, model))
}

eq_fin <- function (x) {
  return(heq(x, S, model))
}



grad <- function (x, S, model) {
  est <- lvm_ml(x, model, jac = FALSE)
  R <- est$p_implied
  lambda <- est$lambda
  BETA <- est$beta
  GAMMA <- est$gamma
  Sigma <- est$sigma_implied

  J_Sigma_theta <- compute_J_Sigma_theta(model$block_sizes, model$lengths_theta, model$mode,model$dag,
                                         R, model$relation_matrix, lambda, BETA, GAMMA)
  # Compute Gradient
  Gradient <- compute_gradient(Sigma, S, J_Sigma_theta)

  return(as.vector(Gradient))
}

grad_heq <- function(x, S, model) {
  est <- lvm_ml(x, model, jac = FALSE)
  lambda <- est$lambda
  S_composites <- est$S_composites


  J_H <- compute_gradient_constraint(lambda, S_composites, model$block_sizes, model$lengths_theta, model$mode)

  return(as.matrix(J_H))
}




system.time({
sol <- csolnp(pars = init, fn =f ,eq_fn =eq_fin, gr = grad, eq_b = rep(0,r), eq_jac = grad_eq, lower = -10*abs(init), upper = 10*abs(init),
              control = list(trace = 0, tol = tol), use_r_version = FALSE)
})


as.vector(grad(init, S, model))

round(result$pars - sol$pars, 3)


system.time({
res0 <- nloptr(x0=init,
        eval_f=f,
        eval_grad_f=grad,
        eval_g_eq = eq_fin,
        eval_jac_g_eq = grad_eq,
        opts = list("algorithm"="NLOPT_LD_SLSQP", 'xtol_rel' = 1e-8))

  })


system.time({
res0bis <- nloptr(x0=init,
        eval_f=F1,
        eval_grad_f=grad,
        eval_g_eq = heq,
        eval_jac_g_eq = grad_eq,
        opts = list("algorithm"="NLOPT_LD_SLSQP", 'xtol_rel' = 1e-8),
        S = S, model = model)

  })





system.time({
res1 <- nloptr(x0=init,
        eval_f=f,
        eval_g_eq = eq_fin,
        opts = list("algorithm"="NLOPT_GN_ISRES", 'xtol_rel' = 1e-8))

})


system.time({
res_slsqp <- nloptr(x0 = init,
                    eval_f = f,
                    eval_grad_f = grad,
                    eval_g_eq = eq_fin,
                    eval_jac_g_eq = grad_eq,
                    opts = list("algorithm" = "NLOPT_LD_SLSQP",
                                "xtol_rel" = 1e-6,
                                "maxeval" = 1000))
})

# On doit définir un sous-algorithme pour l'optimisation interne
local_opts <- list("algorithm" = "NLOPT_LD_LBFGS",
                   "xtol_rel"  = 1.0e-6)
system.time({
res_auglag <- nloptr(x0 = init,
                     eval_f = f,
                     eval_grad_f = grad,
                     eval_g_eq = eq_fin,
                     eval_jac_g_eq = grad_eq,
                     opts = list("algorithm" = "NLOPT_LD_AUGLAG_EQ", # EQ = spécifique pour égalités
                                 "xtol_rel" = 1.0e-6,
                                 "maxeval" = 1000,
                                 "local_opts" = local_opts)) # Ajout du sous-algorithme
})




library(alabama)

# Exécution de l'optimisation
res_alabama <- auglag(
  par = init,               # Point de départ (équivalent à x0)
  fn = f,                   # Fonction objectif (équivalent à eval_f)
  gr = grad,                # Gradient objectif (équivalent à eval_grad_f)
  heq = eq_fin,             # Contraintes d'égalité (équivalent à eval_g_eq)
  heq.jac = grad_eq,        # Jacobien des contraintes (équivalent à eval_jac_g_eq)
  control.outer = list(
    eps = 1e-8,             # Tolérance pour la convergence
    trace = FALSE            # TRUE pour voir la progression s'afficher
  )
)
system.time({
res_constr <- constrOptim.nl(
  par = init,                 # Votre point de départ
  fn = f,                     # Votre fonction objectif
  heq = eq_fin,               # Vos contraintes d'égalité (h(x) = 0)
  control.outer = list(eps = 1e-8, trace = TRUE) # Tolérance et affichage
)
})

round(res0$solution - res1$solution, 3)

round(result$pars - res_alabama$par, 3)
round(result$pars - res_constr$par, 3)

round(res_alabama$par - res_constr$par, 3)
round(res_auglag$solution - res_constr$par, 3)

round(result$pars -res0$solution, 3)


round(result$pars -res_slsqp$solution, 3)


round(res0$solution - res_slsqp$solution, 3)

round(res_auglag$solution - res_slsqp$solution, 3)


round(result$pars -sol$pars, 3)