devtools::load_all()
# remotes::install_github("yrosseel/lavaan")
library(lavaan)
library('Rsolnp')


source("tests/testthat/dataset_configs/load_configs.R")
source("inst/report/report_utils.R")
source("inst/examples/script_gradient/utils.num.R")

dataset_name = 'PoliticalDemocracy'
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
                fun=F1, S = S, model = model,
                control = list(trace = 0, tol = tol))
})
cat("f solnp:", f(result$pars), "\n")

system.time({
    res0 <- nloptr(
      x0=init,
      eval_f=f,
      eval_grad_f=grad,
      lb = rep(-10, length(init)), # Exemple de borne inférieure
       ub = rep(10, length(init)),  # Exemple de borne supérieure
      opts = list("algorithm"="NLOPT_LD_LBFGS",
                  'xtol_rel' = tol)
    )

  })
cat("f nloptr:", f(res0$solution), "\n")

system.time({
res <- nlminb(start = init,
              objective = f,
              gradient = grad,   # Optionnel mais recommandé
              hessian = hess,    # Optionnel mais recommandé
              lower = -Inf,      # Bornes inférieures
              upper = Inf,       # Bornes supérieures
              control = list(x.tol = tol))
  })
cat("f nlminb:", f(res$par), "\n")

system.time({
sol <- csolnp(pars = init, fn =f , gr = grad, lower = -10*abs(init), upper = 10*abs(init),
              control = list(trace = 0, tol = tol), use_r_version = TRUE)
})
cat("f csolnp:", f(sol$pars), "\n")


cat("solnp vs nloptr:", round(result$pars - res0$solution, 3), "\n")
cat("solnp vs csolnp:", round(result$pars - sol$pars, 3), "\n")
cat("nloptr vs csolnp:", round(res0$solution - sol$pars, 3), "\n")
cat("solnp vs nlminb:", round(result$pars - res$par, 3), "\n")


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