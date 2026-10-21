devtools::load_all()
# remotes::install_github("yrosseel/lavaan")
library(lavaan)
library('Rsolnp')


source("tests/testthat/dataset_configs/load_configs.R")
source("inst/report/report_utils.R")
source("inst/examples/script_gradient/utils.num.R")

dataset_name = 'ECSI_full'
dataset <- load_dataset_data(dataset_name)
config <- get_dataset_config(dataset_name, dataset)

svd_results <- fit_and_get_estimates(config, "svd", list(B = 100))
ml_results <- fit_and_get_estimates(config, "ml", list(tol = 1e-8))
system.time({
lavaan_ml <- sem(config$sem,
               data = data.frame(Reduce("cbind", config$data)),
                       estimator = "ML",
                       likelihood = "wishart")

})
lavaan_results <- parameterEstimates(lavaan_ml, standardized = TRUE)


print(merge_estimates_with_lavaan(ml_results$param, lavaan_results))






set.seed(27)
modelsvd <- SemFC$new(data=config$data, relation_matrix = config$relation_matrix, mode=config$mode, estimator = "svd")

modelsvd$fit(infer = F)








init <- modelsvd$get_estimate('theta')

model <- modelsvd$get_model()
dat <- modelsvd$get_data()
S <- dat$cov_S
r = sum(model$mode == 'formative')
tol = 1e-08


cat('f init:', f(init), '\n')
system.time({
result <- solnp(pars = init,
                fun=F1, eqfun=heq1,
                eqB = rep(0,r), S = S, model = model,
                control = list(trace = 0, tol = tol))
})
cat("f solnp:", f(result$pars), "\n")

system.time({
res0 <- nloptr(x0=init,
        eval_f=f,
        eval_grad_f=grad,
        eval_g_eq = h_eq,
        eval_jac_g_eq = grad_h_eq,
        opts = list("algorithm"="NLOPT_LD_SLSQP", 'xtol_rel' = 1e-10))

  })

cat("f nloptr:", f(res0$solution), "\n")
system.time({
sol <- csolnp(pars = init, fn =f ,eq_fn =h_eq, gr = grad, eq_b = rep(0,r), eq_jac = grad_h_eq, lower = -10*abs(init), upper = 10*abs(init),
              control = list(trace = 0, tol = tol), use_r_version = TRUE)
})
cat("f csolnp:", f(sol$pars), "\n")


cat("solnp vs nloptr:", round(result$pars - res0$solution, 3), "\n")
cat("solnp vs csolnp:", round(result$pars - sol$pars, 3), "\n")
cat("nloptr vs csolnp:", round(res0$solution - sol$pars, 3), "\n")


cat("Initial values vs nloptr:", round(init - res0$solution, 3), "\n")



system.time({
grad_analytique <- grad(init)                 # Votre calcul
})

system.time({
grad_numerique  <- numDeriv::grad(f, init)    # Calcul de R par différences finies
})
# Afficher la différence maximale absolue
erreur_max_f <- max(abs(grad_analytique - grad_numerique))
cat("Erreur max pour f :", erreur_max_f, "\n")





# ---------------------------------------------------------
# 2. Vérification des contraintes d'égalité (eq_fin)
# ---------------------------------------------------------
system.time({
jacob_analytique <- grad_h_eq(init)             # Votre calcul
})
system.time({
jacob_numerique  <- numDeriv::jacobian(h_eq, init ) # Calcul de R
})

# Comme vous manipulez des matrices creuses, forcez la conversion
# en matrice classique pour la comparaison
jacob_analytique <- as.matrix(jacob_analytique)

# Afficher la différence maximale absolue
erreur_max_eq <- max(abs(jacob_analytique - jacob_numerique))
cat("Erreur max pour eq_fin :", erreur_max_eq, "\n")






# ---------------------------------------------------------
# 3. check inf
# ---------------------------------------------------------
system.time({
jacob_analytique <- P_ml(init, S, model)             # Votre calcul
})
system.time({
jacob_numerique  <- P_ml_old(init, S, model) # Calcul de R
})

# Comme vous manipulez des matrices creuses, forcez la conversion
# en matrice classique pour la comparaison
jacob_analytique <- as.matrix(jacob_analytique)

# Afficher la différence maximale absolue
erreur_max_eq <- max(abs(jacob_analytique - jacob_numerique))
cat("Erreur max matrice variance covariance :", erreur_max_eq, "\n")