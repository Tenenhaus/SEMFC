

######################################
# Load useful packages and functions #
######################################

source('R/SEMFC/sem_f_c.R')
library(cSEM)


#############################################
########## MONTE-CARLO SIMULATION ###########
#############################################
set.seed(20091979) #my date of birth
n_simu <- 5
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
sigma_hat <- matrix(0, 12, n_simu)
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

      source('inst/model/model_mixed.R')
      # Y <- Y_2
      # X <- X_2

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

      # Likelihood Ratio Test
      # fit.ml0 <- solnp(pars = init_ml_with_S,
      #                 fun=F1, eqfun=heq0,
      #                 eqB = rep(0,5), S = model$S,
      #                 control = list(trace = 0, tol = 1e-8))
      #
      # # ML output
      # x0 <- fit.ml0$pars
      #
      # lrt_pval[b] <- pchisq((N-1)*(
      #    F1(x0, model_ml$S, model_ml$block_sizes, model_ml$mode, model_ml$lengths_parameter, model_ml$which_exo_endo) -
      #      model_ml$parameters$F
      # ), 1, lower.tail = F)

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

      # fit.csem <- csem(.model = sem.model, .data = X)

      fit.csem <- csem(.data = X,
                 .model = sem.model,
                 .approach_weights = "PLS-PM",
                 .PLS_weight_scheme_inner = "factorial",
                 .approach_paths = "2SLS",
                 .instruments = list( eta5 = c("eta1", "eta2", "eta3", "eta4"),
                                      eta6 = c("eta1", "eta2", "eta3", "eta4")),
                 .PLS_ignore_structural_model = TRUE, .tolerance = 1e-8,
                 .disattenuate = TRUE)


      lambda_CSEM <- fit.csem$Estimates$Loading_estimates[fit.csem$Estimates$Loading_estimates!=0]*apply(X, 2, sd)
      omega_CSEM <- fit.csem$Estimates$Weight_estimates
      beta_CSEM <- fit.csem$Estimates$Path_estimates[5:6,5:6]
      gamma_CSEM <- fit.csem$Estimates$Path_estimates[5:6,1:4]





      P_EXO_CSEM <- fit.csem$Estimates$Construct_VCV[1:4,1:4]
      P_ENDO_CSEM <- fit.csem$Estimates$Construct_VCV[5:6,5:6]
      psi_CSEM <- (diag(2)-beta_CSEM)%*%P_ENDO_CSEM%*%t(diag(2)-beta_CSEM) - gamma_CSEM%*%P_EXO_CSEM%*%t(gamma_CSEM)
      residual_variance_CSEM <- list(
        LV5 =
        apply(Y[[5]], 2, var)-(fit.csem$Estimates$Loading_estimates[5, 13:15]*apply(X[, 13:15], 2, sd))^2,
        LV6 =
        apply(Y[[6]], 2, var)-(fit.csem$Estimates$Loading_estimates[6, 16:18]*apply(X[, 16:18], 2, sd))^2
      )



      I_B_1_CSEM <- solve(diag(2)-beta_CSEM)
      r2_1_hat_CSEM <- 1 - I_B_1_CSEM[1, ]%*%psi_CSEM%*%I_B_1_CSEM[1, ]
      r2_2_hat_CSEM <- 1 - I_B_1_CSEM[2, ]%*%psi_CSEM%*%I_B_1_CSEM[2, ]
      R2_CSEM <- c(r2_1_hat_CSEM, r2_2_hat_CSEM)


      R_LVM_11_CSEM <- P_EXO_CSEM
      R_LVM_12_CSEM <- P_EXO_CSEM%*%t(gamma_CSEM)%*%t(I_B_1_CSEM)
      R_LVM_21_CSEM <- I_B_1_CSEM%*%gamma_CSEM%*%t(P_EXO_CSEM)
      R_LVM_22_CSEM <- I_B_1_CSEM%*%(gamma_CSEM%*%P_EXO_CSEM%*%t(gamma_CSEM) + psi_CSEM)%*%t(I_B_1_CSEM)
      R_LVM_CSEM <- rbind(cbind(R_LVM_11_CSEM, R_LVM_12_CSEM),
                         cbind(R_LVM_21_CSEM, R_LVM_22_CSEM))




      parameter_csem <- parameters_svd(lambda = lambda_CSEM,
                            P_EXO = P_EXO_CSEM,
                            G = gamma_CSEM,
                            B = beta_CSEM,
                            P_ENDO = P_ENDO_CSEM,
                            residual_variance = residual_variance_CSEM,
                            S_composites = model$S_composites,
                            mode = model$mode)

      SIGMA_CSEM <- lvm_ml(x = parameter_csem,
                           block_sizes = model$block_sizes,
                           mode = model$mode,
                           lengths_parameter = model$lengths_theta,
                           which_exo_endo = model$which_exo_endo,
                           jac = F,
                           varnames = model$varnames)$SIGMA_IMPLIED



      lambda_hat_csem[, b] <- lambda_CSEM
      omega_hat_csem[, b] <- colSums(omega_CSEM, na.rm = TRUE)[1:12]
      beta_hat_csem[, b] <- beta_CSEM[beta_CSEM!=0]
      gamma_hat_csem[, b] <- gamma_CSEM[gamma_CSEM!=0]
      r2_hat_csem[, b] <- R2_CSEM
      psi_hat_csem[, b] <- psi_CSEM[upper.tri(psi_CSEM, diag = TRUE)]
      var_hat_csem[, b] <- Reduce("c", residual_variance_CSEM)
      rho_hat_csem[, b] <- R_LVM_CSEM[upper.tri(R_LVM_CSEM)]



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










