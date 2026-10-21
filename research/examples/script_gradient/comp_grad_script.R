


library(devtools)
load_all()

set.seed(27)
source('inst/model/model_mixed.R')
Y <- Y_2
X <- X_2

set.seed(27)
model <- SemFC$new(data=Y, relation_matrix = C, mode=mode, estimator = "svd")
system.time({
  model$fit(infer = F)



# model$fit(infer = T, B = 100)
# model$fit(infer = F)




mod = model$get_model()
dat = model$get_data()

block_sizes = c(3,3,3,3,3,3)
lengths_theta = c(18,6,4,2,1,30)
lambda = model$get_estimate('lambda')
Lambda = Matrix::bdiag(lambda)

rel_matrix <- matrix(c(0, 0, 0, 0, 1, 0,
             0, 0, 0, 0, 1, 0,
             0, 0, 0, 0, 0, 1,
             0, 0, 0, 0, 0, 1,
             0, 0, 0, 0, 0, 1,
             0, 0, 0, 0, 1, 0), 6, 6, byrow = TRUE)


BETA = model$get_estimate('beta')
GAMMA = model$get_estimate('gamma')
C = solve(diag(2) - model$get_estimate('beta'))
PHI = model$get_estimate('p_implied')[1:4, 1:4]
PSI = model$get_estimate('psi')
PHI_endo = model$get_estimate('p_implied')[5:6, 5:6]
S_composites = dat$S_diag_composites

h_eq = heq(model$get_estimate('theta'), cov(X), mod)
J_h_eq = compute_gradient_constraint(lambda, S_composites, block_sizes, lengths_theta, mode)

theta_svd = model$get_estimate('theta')




res = compute_J_Sigma_theta(block_sizes, lengths_theta, mode, mod$dag, model$get_estimate('p_implied'),
                                   rel_matrix, model$get_estimate('lambda'), BETA, GAMMA)

H3 = compute_Hessian(model$get_estimate('sigma_implied'), res)
grad_F = as.vector(compute_gradient(model$get_estimate('sigma_implied'), cov(X), res))



Psi_point <- rbind(cbind(H3, t(J_h_eq)),
             cbind(J_h_eq, matrix(0, 4, 4)))

Psi = c(grad_F, h_eq)


theta_os = theta_svd - solve(Psi_point, Psi)[1:length(theta_svd)]







})

loadings_os_form <- get_loadings(theta_os, block_sizes)[mode == "formative"]


BDIAG <- get_bdiag(theta_os,
                   mode = mode,
                   block_sizes = block_sizes,
                   initial_start_index_cov = cumsum(c(1, head(mod$lengths_theta, -1)))[6])

S_composites_os <- BDIAG[mode == 'formative']



loadings_final_form <- Map(function(x, M) {
  # 1. Calcul du dénominateur (la norme de Mahalanobis)
  # as.numeric() garantit que l'on obtient un scalaire et non une matrice 1x1
  norm_factor <- sqrt(as.numeric(t(x) %*% solve(M, x)))

  # 2. Rétraction : division du vecteur par le scalaire
  return(x / norm_factor)

}, loadings_os_form, S_composites_os)




theta_final <- theta_os
theta_final[1:12] <- unlist(loadings_final_form)


system.time({
modelml <- SemFC$new(data=Y, relation_matrix = C, mode=mode, estimator = "ml")
modelml$fit(infer = T)
})
theta_ml = modelml$get_estimate('theta')

round(theta_os - theta_ml, 3)
round(theta_os - theta_svd, 3)



set.seed(27)
source('inst/model/model_reflective.R')
Y <- Y_2
X <- X_2