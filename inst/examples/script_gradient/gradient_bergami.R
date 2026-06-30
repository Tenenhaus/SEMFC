
library(devtools)
load_all()

set.seed(27)
source('inst/model/model_bergami.R')


set.seed(27)
model <- SemFC$new(data=A_bergami, relation_matrix = C_bergami, mode=mode_bergami, estimator = "svd")


model$fit(infer = F)




mod = model$get_model()
mode = mod$mode

block_sizes = mod$block_sizes
lengths_theta = mod$lengths_theta
lambda = model$get_estimate('lambda')
Lambda = Matrix::bdiag(lambda)

rel_matrix <- mod$relation_matrix


BETA = model$get_estimate('beta')
GAMMA = model$get_estimate('gamma')
C = solve(diag(ncol(BETA)) - BETA)
Phi_exo = model$get_estimate('p_exo')
PSI = model$get_estimate('psi')
Phi_endo = model$get_estimate('p_endo')
R = model$get_estimate('p_implied')





system.time({
  I = information_matrix(model$get_estimate('theta'), mod)
})



system.time({
res = compute_J_Sigma_theta(block_sizes, lengths_theta, mode, mod$dag, model$get_estimate('p_implied'),
                                   rel_matrix, model$get_estimate('lambda'), BETA, GAMMA)

H3 = compute_Hessian(model$get_estimate('sigma_implied'), res)

I2 = H3 * 0.5
})
#


sum(I - I2)