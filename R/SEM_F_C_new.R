##################################################################
# Reproducible results for Tenenhaus, Tenenhaus & Dijkstra paper #
#              submitted to ADAC in July 15th, 2024              #
##################################################################

#############################################
# Remove all objects from the R environment #
#############################################
rm(list = ls())

######################################
# Load useful packages and functions #
######################################

source('R/SEMFC/sem_f_c.R')






#############################################
########## MONTE-CARLO SIMULATION ###########
#############################################
set.seed(20091979) #my date of birth
n_simu <- 100
N <- 300
sol_svd <- matrix(0, 61, n_simu)
sol_ml <- matrix(0, 61, n_simu)


#eigen
lambda_hat_svd <- matrix(0, 18, n_simu)
omega_hat_svd <- matrix(0, 12, n_simu)
rho_hat_svd <- matrix(0, 15, n_simu)
beta_hat_svd <- matrix(0, 2, n_simu)
gamma_hat_svd <- matrix(0, 4, n_simu)
r2_hat_svd <- matrix(0, 2, n_simu)
psi_hat_svd <- matrix(0, 3, n_simu)
var_hat_svd <- matrix(0, 6, n_simu)
f_svd <- rep(0, n_simu)
param_svd <- matrix(0, 61, n_simu)
z <- rep(NA, n_simu)
svd_pval <- rep(NA, n_simu)

#Maximum likelihood
lambda_hat_ml <- matrix(0, 18, n_simu)
omega_hat_ml <- matrix(0, 12, n_simu)
rho_hat_ml <- matrix(0, 15, n_simu)
beta_hat_ml <- matrix(0, 2, n_simu)
gamma_hat_ml <- matrix(0, 4, n_simu)
r2_hat_ml <- matrix(0, 2, n_simu)
psi_hat_ml <- matrix(0, 3, n_simu)
sigma_hat <- matrix(0, 12, n_simu)
var_hat_ml <- matrix(0, 6, n_simu)
std_err_ml <- matrix(NA, 61, n_simu)
f_ml <- rep(0, n_simu)
param_ml <- matrix(0, 61, n_simu)
lrt_pval <- rep(NA, n_simu)
vcov_pval <- rep(NA, n_simu)

# Goodness of fit
gof <- matrix(NA, 2, n_simu)


for (b in seq_len(n_simu)){
  if (b%%10==0) print(b)
  try(
    {

      source('data/data_generated_mixed.R')
      Y <- Y_2

      model <- SemFC$new(data=Y, relation_matrix = C, mode=mode, scale=F, bias=F)
      model$fit_svd()
      model$svd_infer()


      z[b] <- diff(model$infer_estimate$gamma[2:3,1])/sd(diff(t(model$infer_estimate$out[[4]][, 2:3])))
      svd_pval[b] <- 2*pnorm(abs(z[b]), lower.tail = FALSE)


      #################
      # Sigma_implied #
      #################

      lambda_SVD <- model$parameters$lambda
      omega_SVD <- model$parameters$omega
      beta_SVD <- model$parameters$beta
      gamma_SVD <- model$parameters$gamma
      R2_SVD <- model$parameters$R2
      psi_SVD <- model$parameters$psi
      residual_variance_SVD <- model$parameters$residual_variance

      R_LVM_SVD <- model$parameters$P_IMPLIED
      SIGMA_SVD <- model$parameters$SIGMA_IMPLIED


      lambda_hat_svd[, b] <- Reduce("c", lambda_SVD)
      omega_hat_svd[, b] <- Reduce("c", omega_SVD)
      rho_hat_svd[, b] <- R_LVM_SVD[as.vector(upper.tri(R_LVM_SVD))]
      beta_hat_svd[, b] <- beta_SVD[beta_SVD!=0]
      gamma_hat_svd[, b] <- gamma_SVD[gamma_SVD!=0]
      r2_hat_svd[, b] <- R2_SVD
      psi_hat_svd[, b] <- psi_SVD[upper.tri(psi_SVD, diag = TRUE)]
      var_hat_svd[, b] <- Reduce("c", residual_variance_SVD)

      init_ml_with_S <- model$parameters$theta

      sol_svd[, b] <- init_ml_with_S

      # ML output
      model_ml <- SemFC$new(data=Y, relation_matrix = C, mode=mode, scale=F, bias=F)
      model_ml$fit_ml(initialisation_svd  = TRUE)
      model_ml$ml_infer()

      x <- model_ml$parameters$theta
      VCOV <- model_ml$VCOV
      SD <- model_ml$SD

      z <- abs(diff(x[26:27]))/sqrt(VCOV[26, 26] + VCOV[27, 27] - 2*VCOV[26, 27])
      vcov_pval[b] <- 2*pnorm(z, lower.tail = F)

      P_EXO_ML <- model_ml$parameters$P_EXO
      P_ENDO_ML <- model_ml$parameters$P_ENDO
      # path coefficients
      gamma_ML <- model_ml$parameters$gamma
      # path coefficients
      beta_ML <- model_ml$parameters$beta
      R_LVM_ML <- model_ml$parameters$R_LVM
      psi_ML <-  model_ml$parameters$psi
      R2_ML <- model_ml$parameters$R2
      SIGMA_ML <- model_ml$parameters$SIGMA_IMPLIED
      omega_ML <- model_ml$parameters$omega
      lambda_ML <- model_ml$parameters$lambda
      residual_variance_ML <- model_ml$parameters$residual_variance


      lambda_hat_ml[, b] <- Reduce('c', lambda_ML)
      omega_hat_ml[, b] <-  Reduce('c', omega_ML)
      rho_hat_ml[, b] <- R_LVM_ML[upper.tri(R_LVM_ML)]
      beta_hat_ml[, b] <- beta_ML[beta_ML!=0]
      gamma_hat_ml[, b] <- gamma_ML[gamma_ML!=0]
      r2_hat_ml[, b] <- R2_ML
      psi_hat_ml[, b] <- psi_ML[upper.tri(psi_ML, diag = TRUE)]
      var_hat_ml[, b] <- Reduce('c', residual_variance_ML)
      std_err_ml[, b] <- SD


      # SIGMA
      S_composites_empirical <-model$S_composites
      S_composites_ML <- model_ml$parameters$S_composites
      S_composites_true <- list(SIGMA11, SIGMA22,SIGMA33,SIGMA44)

      dls_empirical_true <- mapply(d_LS, S_composites_empirical,S_composites_true, SIMPLIFY = T)
      dls_ml_true <- mapply(d_LS, S_composites_ML,S_composites_true, SIMPLIFY = T)
      dls_ml_empirical <- mapply(d_LS, S_composites_ML,S_composites_empirical, SIMPLIFY = T)

      sigma_hat[, b] <- c(dls_empirical_true, dls_ml_true, dls_ml_empirical)

      #goodness of fit
      gof[1, b] <- d_LS(SIGMA_SVD, SIGMA)
      gof[2, b] <- d_LS(SIGMA_ML, SIGMA)


      # param eigen
      param_svd[, b] <- init_ml_with_S

      # param ml
      param_ml[, b] <- x


      # F_ml
      f_ml[b] <- model_ml$parameters$F

      # F_svd
      f_svd[b] <- model$parameters$F


    }, silent = TRUE
  )

}



