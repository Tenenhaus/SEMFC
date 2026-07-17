devtools::load_all()
set.seed(27)



source("tests/testthat/dataset_configs/load_configs.R")
source("inst/report/report_utils.R")
source("inst/examples/script_gradient/utils.num.R")

# choisir un dataset formatif
dataset_name = 'ECSI_full'



dataset <- load_dataset_data(dataset_name)
config <- get_dataset_config(dataset_name, dataset)


modelsvd <- SemFC$new(data=config$data, relation_matrix = config$relation_matrix, mode=config$mode, estimator = "svd")

modelsvd$fit(infer = F)








init <- modelsvd$get_estimate('theta')
model <- modelsvd$get_model()
dat <- modelsvd$get_data()
S <- dat$cov_S
r = sum(model$mode == 'formative')
tol = 1e-08





BETA = modelsvd$get_estimate('beta')
GAMMA = modelsvd$get_estimate('gamma')
S_composites = dat$S_diag_composites
lambda = modelsvd$get_estimate('lambda')
block_sizes = model$block_sizes
lengths_theta = model$lengths_theta
S = dat$cov_S

h_eq = heq(modelsvd$get_estimate('theta'), cov(X), model)
J_h_eq = compute_gradient_constraint(lambda, S_composites, block_sizes, lengths_theta, config$mode)

theta_svd = modelsvd$get_estimate('theta')




res = compute_J_Sigma_theta(block_sizes, lengths_theta, config$mode, model$dag, modelsvd$get_estimate('p_implied'),
                                   config$relation_matrix, lambda, BETA, GAMMA)

H3 = compute_Hessian(modelsvd$get_estimate('sigma_implied'), res)
grad_F = as.vector(compute_gradient(modelsvd$get_estimate('sigma_implied'), S, res))



Psi_point <- rbind(cbind(H3, t(J_h_eq)),
             cbind(J_h_eq, matrix(0, r, r)))

Psi = c(grad_F, h_eq)


theta_os = theta_svd - solve(Psi_point, Psi)[1:length(theta_svd)]


#normalisation des loadings sujets aux contraintes des blocs formatifs
loadings_os_form <- get_loadings(theta_os, block_sizes)[config$mode == "formative"]


BDIAG <- get_bdiag(theta_os,
                   mode = config$mode,
                   block_sizes = block_sizes,
                   initial_start_index_cov = cumsum(c(1, head(model$lengths_theta, -1)))[6])

S_composites_os <- BDIAG[config$mode == 'formative']



loadings_final_form <- Map(function(x, M) {

  norm_factor <- sqrt(as.numeric(t(x) %*% solve(M, x)))


  return(x / norm_factor)

}, loadings_os_form, S_composites_os)





theta_final <- theta_os
theta_final[1:length(unlist(loadings_final_form))] <- unlist(loadings_final_form)


cat("f init:", f(init), '\n')
cat('f theta_os:', f(theta_os), '\n')
cat('f theta_final:', f(theta_final), '\n')


modelml <- SemFC$new(data=config$data, relation_matrix = config$relation_matrix, mode=config$mode, estimator = "ml")

modelml$fit(infer = F)


cat('f theta_ml:', f(modelml$get_estimate('theta')), '\n')




