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
library(cSEM)


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


# csem

lambda_hat_csem <- matrix(0, 18, n_simu)
omega_hat_csem <- matrix(0, 12, n_simu)
rho_hat_csem <- matrix(0, 15, n_simu)
beta_hat_csem <- matrix(0, 2, n_simu)
gamma_hat_csem <- matrix(0, 4, n_simu)
r2_hat_csem <- matrix(0, 2, n_simu)
psi_hat_csem <- matrix(0, 3, n_simu)
var_hat_csem <- matrix(0, 6, n_simu)
f_csem <- rep(0, n_simu)
param_csem <- matrix(0, 61, n_simu)
z_csem <- rep(NA, n_simu)
csem_pval <- rep(NA, n_simu)



#Maximum likelihood
lambda_hat_ml <- matrix(0, 18, n_simu)
omega_hat_ml <- matrix(0, 12, n_simu)
rho_hat_ml <- matrix(0, 15, n_simu)
beta_hat_ml <- matrix(0, 2, n_simu)
gamma_hat_ml <- matrix(0, 4, n_simu)
r2_hat_ml <- matrix(0, 2, n_simu)
psi_hat_ml <- matrix(0, 3, n_simu)
sigma_hat <- matrix(0, 20, n_simu)
var_hat_ml <- matrix(0, 6, n_simu)
std_err_ml <- matrix(NA, 61, n_simu)
f_ml <- rep(0, n_simu)
param_ml <- matrix(0, 61, n_simu)
lrt_pval <- rep(NA, n_simu)
vcov_pval <- rep(NA, n_simu)

# Goodness of fit
gof <- matrix(NA, 3, n_simu)


for (b in seq_len(n_simu)){
  if (b%%10==0) print(b)
  try(
    {

      source('data/data_generated_mixed.R')
      Y <- Y_2
      X <- X_2

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


      # csem

      fit.csem <- csem(.model = sem.model, .data = X)


      lambda_CSEM <- fit.csem$Estimates$Loading_estimates
      omega_CSEM <- fit.csem$Estimates$Weight_estimates
      beta_CSEM <- fit.csem$Estimates$Path_estimates[5:6,5:6]
      gamma_CSEM <- fit.csem$Estimates$Path_estimates[5:6,1:4]
      R2_CSEM <- fit.csem$Estimates$R2

      P_EXO_CSEM <- fit.csem$Estimates$Construct_VCV[1:4,1:4]
      P_ENDO_CSEM <- fit.csem$Estimates$Construct_VCV[5:6,5:6]
      psi_CSEM <- (diag(2)-beta_CSEM)%*%P_ENDO_CSEM%*%t(diag(2)-beta_CSEM) - gamma_CSEM%*%P_EXO_CSEM%*%t(gamma_CSEM)
      residual_variance_CSEM <- diag(fit.csem$Estimates$Residual_correlation)[13:18]

      R_LVM_CSEM  <- fit.csem$Estimates$Construct_VCV
      SIGMA_CSEM <- fit.csem$Estimates$Indicator_VCV

      SIGMA11_CSEM <- SIGMA_CSEM[1:3, 1:3]
      SIGMA22_CSEM <- SIGMA_CSEM[4:6, 4:6]
      SIGMA33_CSEM <- SIGMA_CSEM[7:9, 7:9]
      SIGMA44_CSEM <- SIGMA_CSEM[10:12, 10:12]


      parameter_csem <- Reduce("c",
                                c(Reduce("c", colSums(lambda_CSEM, na.rm = TRUE)),
                                  P_EXO_CSEM[upper.tri(P_EXO_CSEM)],
                                  apply(gamma_CSEM, 1, function(row) row[row != 0]),
                                  apply(beta_CSEM, 1, function(row) row[row != 0]),
                                  P_ENDO_CSEM[upper.tri(P_ENDO_CSEM)],
                                  SIGMA11_CSEM[upper.tri(SIGMA11_CSEM, diag = TRUE)],
                                  SIGMA22_CSEM[upper.tri(SIGMA22_CSEM, diag = TRUE)],
                                  SIGMA33_CSEM[upper.tri(SIGMA33_CSEM, diag = TRUE)],
                                  SIGMA44_CSEM[upper.tri(SIGMA44_CSEM, diag = TRUE)],
                                  residual_variance_CSEM
                                )
    )



      lambda_hat_csem[, b] <- colSums(lambda_CSEM, na.rm = TRUE)
      omega_hat_csem[, b] <- colSums(omega_CSEM, na.rm = TRUE)[1:12]
      beta_hat_csem[, b] <- beta_CSEM[beta_CSEM!=0]
      gamma_hat_csem[, b] <- gamma_CSEM[gamma_CSEM!=0]
      r2_hat_csem[, b] <- R2_CSEM
      psi_hat_csem[, b] <- psi_CSEM[upper.tri(psi_CSEM, diag = TRUE)]
      var_hat_csem[, b] <- residual_variance_CSEM



      # SIGMA
      S_composites_empirical <-model$S_composites
      S_composites_ML <- model_ml$parameters$S_composites
      S_composites_CSEM <- list(SIGMA11_CSEM,
                                SIGMA22_CSEM,
                                SIGMA33_CSEM,
                                SIGMA44_CSEM)
      S_composites_true <- list(SIGMA11, SIGMA22,SIGMA33,SIGMA44)

      dls_empirical_true <- mapply(d_LS, S_composites_empirical,S_composites_true, SIMPLIFY = T)

      dls_csem_true <- mapply(d_LS, S_composites_CSEM,S_composites_true, SIMPLIFY = T)
      dls_ml_true <- mapply(d_LS, S_composites_ML,S_composites_true, SIMPLIFY = T)
      dls_csem_empirical <- mapply(d_LS, S_composites_CSEM,S_composites_empirical, SIMPLIFY = T)
      dls_ml_empirical <- mapply(d_LS, S_composites_ML,S_composites_empirical, SIMPLIFY = T)

      sigma_hat[, b] <- c(dls_empirical_true, dls_ml_true, dls_ml_empirical, dls_csem_true, dls_csem_empirical)

      #goodness of fit
      gof[1, b] <- d_LS(SIGMA_SVD, SIGMA)
      gof[2, b] <- d_LS(SIGMA_ML, SIGMA)
      gof[3, b] <- d_LS(SIGMA_CSEM, SIGMA)


      # param eigen
      param_svd[, b] <- init_ml_with_S

      # param ml
      param_ml[, b] <- x

      # param csem
      param_csem[, b] <- parameter_csem


      # F_ml
      f_ml[b] <- model_ml$parameters$F

      # F_svd
      f_svd[b] <- model$parameters$F

      # F_csem
      f_csem[b] <-  F1(parameter_csem, model$cov_S, model$block_sizes,
                       model$mode, model$lengths_theta, model$which_exo_endo)


    }, silent = TRUE
  )

}

source('R/table_parameters.R')









