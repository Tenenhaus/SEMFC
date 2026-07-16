


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
# system.time({
#   I = information_matrix(model$get_estimate('theta'), mod)
# })



#
# # system.time({
#   J_Sigma_Lambda = jac_Sigma_Lambda(Lambda, model$get_estimate('p_implied'))
#
#   J_Lambda = calcul_J_Lambda_theta(block_sizes, lengths_theta)
#
#
#   J_Sigma_R = jac_Sigma_R(Lambda)
#
#
#   J_Theta = jac_Theta_theta(lambda, block_sizes, lengths_theta, mode)
#
#
#
#
#
#
#
#   P_endo <- P_exo_endo(m_total = 6, m_exo = 4)$P_endo
#   P_exo <- P_exo_endo(m_total = 6, m_exo = 4)$P_exo
#   M_gamma <- M_beta_gamma(rel_matrix)$M_gamma
#   M_beta <- M_beta_gamma(rel_matrix)$M_beta
#
#   J_rhoR_endo = calcul_J_rhoR_endo(P_endo, M_endo = correlation_duplication_matrix(2), L_bar = correlation_elimination_matrix(6))
#
#
#
#
#   J_rhoR_exo = calcul_J_rhoR_exo(P_endo, P_exo, C, GAMMA, M_exo = correlation_duplication_matrix(4), L_bar = correlation_elimination_matrix(6), K_m = commutation_matrix(6))
#
#
#   J_rhoR_gamma <- calcul_J_rhoR_Gamma(P_endo, P_exo, C, PHI, M_gamma)
#
#   J_rhoR_B <- calcul_J_rhoR_B(P_endo, P_exo, C, GAMMA, PHI, M_beta)
#
#
#   J_rhoR = calcul_J_rhoR_theta(J_rhoR_exo, J_rhoR_gamma, J_rhoR_B, J_rhoR_endo, n_lambda = lengths_theta[1], n_Theta = tail(lengths_theta, 1))
#
#
#
#
#
#
#   # Jac_Sigma_theta = J_Sigma_Lambda %*% J_Lambda + J_Sigma_R %*%  J_rhoR + J_Theta
#   Jac_Sigma_theta = jac_vech_Sigma(J_Sigma_Lambda, J_Lambda, J_Sigma_R, J_rhoR, J_Theta)
#   H3 = calcul_Hessian(model$get_estimate('sigma_implied'), duplication_matrix(18), Jac_Sigma_theta)
#   I2 = H3 * 0.5
#
# # })
#
#
# sum(I - I2)

# system.time({
# res = compute_J_Sigma_theta(block_sizes, lengths_theta, mode, mod$dag, model$get_estimate('p_implied'),
#                                    rel_matrix, model$get_estimate('lambda'), BETA, GAMMA)
#
# H3 = compute_Hessian(model$get_estimate('sigma_implied'), res)
#
# I2 = H3 * 0.5
# })
# #


# sum(I - I2)



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