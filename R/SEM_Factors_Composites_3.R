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

# packages
library(readxl)
library(roxygen2)
library(devtools)
library(Matrix)
library(knitr)
library(pheatmap)
library(mvtnorm)
library(MASS)
library(lavaan)
library(Rsolnp)
# load_all("RGCCA")
source("RGCCA/R/svdSEM.R")
source("RGCCA/R/scale2.R")
source("RGCCA/R/correction.R")
source("RGCCA/R/cov2.R")
source("RGCCA/R/lvm.R")
source("RGCCA/R/ind_exo_endo.R")
source("RGCCA/R/d_LS.R")
source("RGCCA/R/improper.R")
source("RGCCA/R/svdSEM_infer.R")
source("RGCCA/R/scaleDataSet.R")
source("RGCCA/R/svdSEM_gof.R")
source("RGCCA/R/model_sem.R")
source("RGCCA/R/sem_f_c.R")


#functions
source("SEMFC/data/data_simulation.R")
source("SEMFC/ml_sem/F1.R")
source("../functions/h_theta.R")
source("functions/s_implied.R")
source("SEMFC/ml_sem/h_constraints.R")
source("SEMFC/svd_sem/parameters_svd.R")
source("functions/mlSEM.R")
source("functions/rml_sem.R")

#############################################
########## MONTE-CARLO SIMULATION ###########
#############################################
set.seed(20091979) #my date of birth
n_simu = 1000
N = 300
sol_svd = matrix(0, 61, n_simu)
sol_ml = matrix(0, 61, n_simu)

#eigen
lambda_hat_svd = matrix(0, 18, n_simu)
omega_hat_svd = matrix(0, 12, n_simu)
rho_hat_svd = matrix(0, 15, n_simu)
beta_hat_svd = matrix(0, 2, n_simu)
gamma_hat_svd = matrix(0, 4, n_simu)
r2_hat_svd = matrix(0, 2, n_simu)
psi_hat_svd = matrix(0, 3, n_simu)
var_hat_svd = matrix(0, 6, n_simu)
f_svd = rep(0, n_simu)
param_svd = matrix(0, 61, n_simu)
z = rep(NA, n_simu)
svd_pval = rep(NA, n_simu)

#Maximum likelihood
lambda_hat_ml = matrix(0, 18, n_simu)
omega_hat_ml = matrix(0, 12, n_simu)
rho_hat_ml = matrix(0, 15, n_simu)
beta_hat_ml = matrix(0, 2, n_simu)
gamma_hat_ml = matrix(0, 4, n_simu)
r2_hat_ml = matrix(0, 2, n_simu)
psi_hat_ml = matrix(0, 3, n_simu)
sigma_hat = matrix(0, 12, n_simu)
var_hat_ml = matrix(0, 6, n_simu)
std_err_ml = matrix(NA, 61, n_simu)
f_ml = rep(0, n_simu)
param_ml = matrix(0, 61, n_simu)
lrt_pval = rep(NA, n_simu)
vcov_pval = rep(NA, n_simu)

# Goodness of fit
gof = matrix(NA, 2, n_simu)

for (b in seq_len(n_simu)){
  if (b%%10==0) print(b)
  try(
    {
      X = mvrnorm(N, rep(0, 18), SIGMA, empirical = FALSE)
      colnames(X) = paste("X", rep(1:6, each = 3), rep(1:3, 6), sep ="")

      
      Y = list(X1 = X[, 1:3], X2 = X[, 4:6], X3 = X[, 7:9],
               X4 = X[, 10:12], X5 = X[, 13:15], X6 = X[, 16:18])


      
      C = matrix(c(0, 0, 0, 0, 1, 0,
                   0, 0, 0, 0, 1, 0,
                   0, 0, 0, 0, 0, 1,
                   0, 0, 0, 0, 0, 1,
                   0, 0, 0, 0, 0, 1,
                   0, 0, 0, 0, 1, 0), 6, 6, byrow = TRUE)
      
      colnames(C) = rownames(C) = names(Y)

      mode <- c(rep("formative", 4), rep("reflective", 2))

      model <- SEM_F_C$new(data=Y, relation_matrix = C, mode=mode, scale=F, bias=F)
      model$fit_svd()


      model$svd_infer()
      boot_out <- model$boot_svd


      z[b] = diff(boot_out$gamma[2:3,1])/sd(diff(t(boot_out$out[[4]][, 2:3])))
      svd_pval[b] = 2*pnorm(abs(z[b]), lower.tail = FALSE)
      
      #################
      # Sigma_implied #
      #################

      lambda_SVD <- model$svd_parameters$lambda
      omega_SVD <- model$svd_parameters$omega
      beta_SVD <- model$svd_parameters$beta
      gamma_SVD <- model$svd_parameters$gamma
      R2_SVD <- model$svd_parameters$R2
      psi_SVD <- model$svd_parameters$psi
      residual_variance_SVD <- model$svd_parameters$residual_variance

      R_LVM_SVD <- model$svd_parameters$P_IMPLIED
      SIGMA_SVD <- model$svd_parameters$SIGMA_IMPLIED

      
      
      
      lambda_hat_svd[, b] <- Reduce("c", lambda_SVD)
      omega_hat_svd[, b] <- Reduce("c", omega_SVD)
      rho_hat_svd[, b] <- R_LVM_SVD[as.vector(upper.tri(R_LVM_SVD))]
      beta_hat_svd[, b] <- beta_SVD[beta_SVD!=0]
      gamma_hat_svd[, b] <- gamma_SVD[gamma_SVD!=0]
      r2_hat_svd[, b] <- R2_SVD
      psi_hat_svd[, b] <- psi_SVD[upper.tri(psi_SVD, diag = TRUE)]
      var_hat_svd[, b] <- Reduce("c", residual_variance_SVD[mode == 'reflective' ])

      init_ml_with_S <- model$svd_parameters$theta
      
      sol_svd[, b] <- init_ml_with_S


      model$fit_ml(initialisation_svd = TRUE)
      
      # ML output
      x <- model$ml_parameters$theta
      
            
      # # Likelihood Ratio Test
      # fit.ml0 = solnp(pars = init_ml_with_S,
      #                 fun=F1, eqfun=heq0,
      #                 eqB = rep(0,5), S = S,
      #                 control = list(trace = 0, tol = 1e-8))
      #
      # # ML output
      # x0 = fit.ml0$pars
      #
      # lrt_pval[b] = pchisq((N-1)*(F1(x0, S)-F1(x, S)),
      #                       1,
      #                       lower.tail = F)
      #
      #
      # asymptotic test based on the CAN properties of ML 
      # -> asymptotic covariance matrix of theta
      
      # Jacobian and hessian comutation

      model$ml_infer()

      VCOV <- model$VCOV
      SD <- model$SD
      
      z <- abs(diff(x[26:27]))/sqrt(VCOV[26, 26] + VCOV[27, 27] - 2*VCOV[26, 27])
      vcov_pval[b] = 2*pnorm(z, lower.tail = F)

      parameter_implied = s_implied_bis(x, block_sizes, mode, lengths_parameter, which_exo_endo,jac=FALSE)
      P_EXO_ML <- model$ml_parameters$P_EXO
      P_ENDO_ML <- parameter_implied$P_ENDO
      # path coefficients
      gamma_ML <- parameter_implied$GAMMA
      # path coefficients
      beta_ML <- parameter_implied$BETA
      R_LVM_ML <- parameter_implied$R_LVM
      psi_ML <-  parameter_implied$PSI
      R2_ML <- parameter_implied$R2
      SIGMA_ML <- parameter_implied$implied_S
      omega_ML <- parameter_implied$omega
      lambda_ML <- parameter_implied$lambda
      residual_variance_ML <- parameter_implied$residual_variance
      
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

      S_composites_ML <- model$S_composites
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
      f_ml[b] <- model$ml_parameters$F
      
      # F_svd
      f_svd[b] <- model$svd_parameters$F
    }, silent = TRUE
  )
}


##################################################
############ Table IMPROPER SOLUTIONS ############
##################################################

N = c(seq(20, 100, by = 10), 200, 300, 400, 500, 600, 700)
n_simu =  1000
n_improper = 9
n_improper_ml = 7
improper_sol = array(NA, dim = c(n_simu, n_improper, length(N)))
improper_sol_ml = array(NA, dim = c(n_simu, n_improper_ml, length(N)))
Table_improper = matrix(0, 9, length(N))
Table_improper_ml = matrix(0, 7, length(N))


set.seed(27) #my favorite number
for (n in 1:length(N)){
  for(s in 1:n_simu){
    print(paste(N[n], " & ", s))
    X = mvrnorm(N[n], rep(0, 18), SIGMA, empirical = FALSE)
    colnames(X) = paste("X", rep(1:6, each = 3), rep(1:3, 6), sep ="")
    S = cov(X)
    
    Y = list(X1 = X[, 1:3], X2 = X[, 4:6], X3 = X[, 7:9],
             X4 = X[, 10:12], X5 = X[, 13:15], X6 = X[, 16:18])
    
    C = matrix(c(0, 0, 0, 0, 1, 0,
                 0, 0, 0, 0, 1, 0,
                 0, 0, 0, 0, 0, 1,
                 0, 0, 0, 0, 0, 1,
                 0, 0, 0, 0, 0, 1,
                 0, 0, 0, 0, 1, 0), 6, 6, byrow = TRUE)
    
    colnames(C) = rownames(C) = names(Y)
    
    fit.svd = svdSEM(Y, C, scale = FALSE, 
                    mode = c(rep("formative", 4), 
                             rep("reflective", 2)), 
                    bias = FALSE) 
    
    improper_sol[s, , n] = improper(fit.svd)
    
    dimnames(improper_sol) = list(1:n_simu, c("RELIABILITY_COEF",
                                              "RHO_JH", 
                                              "P_TILDE", 
                                              "P_IMPLIED", 
                                              "THETA_JH", 
                                              "STD_LAMBDA", 
                                              "SIGMA_IMPLIED",
                                              "R2",
                                              "PSI"), N)
    
  try({  
    #################################
    # Restricted Maximum Likelihood #
    #################################
    
    init_ml = as.vector(c(Reduce("c", fit.svd$lambda),
                          fit.svd$P_IMPLIED[1:4, 1:4][upper.tri(fit.svd$P_IMPLIED[1:4, 1:4])],
                          fit.svd$gamma[1, 1:2], fit.svd$gamma[2, 3:4],
                          fit.svd$beta[1, 2], fit.svd$beta[2, 1],
                          fit.svd$P_IMPLIED[5,6],
                          S[1:3, 1:3][upper.tri(S[1:3, 1:3], diag = TRUE)],
                          S[4:6, 4:6][upper.tri(S[4:6, 4:6], diag = TRUE)],
                          S[7:9, 7:9][upper.tri(S[7:9, 7:9], diag = TRUE)],
                          S[10:12, 10:12][upper.tri(S[10:12, 10:12], diag = TRUE)],
                          apply(X[, 13:15], 2, var)- fit.svd$lambda[[5]]^2,
                          apply(X[, 16:18], 2, var)- fit.svd$lambda[[6]]^2
    ))
    
    
    fit.ml = solnp(pars = init_ml,
                   fun=F1, eqfun=heq1,
                   eqB = rep(0,4), S = S,
                   control = list(trace = 0, tol = 1e-4))
    
    # ML output
    x = fit.ml$pars
    
    
    
    P_EXO = matrix(c(1, x[19], x[20], x[22],
                     x[19], 1, x[21], x[23],
                     x[20], x[21], 1, x[24],
                     x[22], x[23], x[24], 1), 4, 4)
    
    P_ENDO = matrix(c(1, x[31],
                      x[31], 1), 2, 2)
    
    # path coefficients
    G = matrix(c(x[25], x[26], 0, 0,
                 0, 0, x[27], x[28]), 2, 4, byrow = TRUE)
    
    # path coefficients
    B = matrix(c(0, x[29],
                 x[30], 0), 2, 2, byrow = TRUE)
    
    R_LVM_ML = rbind(cbind(P_EXO, P_EXO%*%t(G)%*%t(solve(diag(2) - B))),
                     cbind(solve(diag(2) - B)%*%G%*%P_EXO, P_ENDO))
    
    I_B_ML = diag(2) - B
    
    PSI_ML =  I_B_ML%*%P_ENDO%*%t(I_B_ML) - G%*%P_EXO%*%t(G)
    
    r2_1_ML = 1 - PSI_ML[1, 1]
    r2_2_ML = 1 - PSI_ML[2, 2 ]
    R2_ML = c(r2_1_ML, r2_2_ML)
    
    #Omega
    #cov between MVs
    S1 = matrix(c(x[32], x[33], x[35],
                  x[33], x[34], x[36],
                  x[35], x[36], x[37]), 3, 3)
    
    #cov between MVs
    
    S2 = matrix(c(x[38], x[39], x[41],
                  x[39], x[40], x[42],
                  x[41], x[42], x[43]), 3, 3)
    
    #cov between MVs
    S3 = matrix(c(x[44], x[45], x[47],
                  x[45], x[46], x[48],
                  x[47], x[48], x[49]), 3, 3)
    #cov between MVs
    S4 = matrix(c(x[50], x[51], x[53],
                  x[51], x[52], x[54],
                  x[53], x[54], x[55]), 3, 3)
    
    LAMBDA_ML = bdiag(split(x[1:18], f = rep(1:6, each = 3)))
    STD_LAMBDA_ML = LAMBDA_ML/apply(X, 2, var)
    RESID_VAR_ML = apply(X, 2, var) - x[1:18]^2
    SIGMA_ML = LAMBDA_ML%*%R_LVM_ML%*%t(LAMBDA_ML) + diag(RESID_VAR_ML)
    
    SIGMA_ML[1:3, 1:3] = S1      #S[1:3, 1:3]
    SIGMA_ML[4:6, 4:6] = S2      #S[4:6, 4:6]
    SIGMA_ML[7:9, 7:9] = S3      #S[7:9, 7:9]
    SIGMA_ML[10:12, 10:12] = S4  #S[10:12, 10:12]
    
    
    
    
    # Improper solution 1 : number of out of range rho_jk 
    # Expected value : -1 <=rho_jk <=1
    
    improper_sol_ml[s, 1, n] = any(abs(R_LVM_ML[upper.tri(R_LVM_ML)])>1)
    
    # Improper solution 2 : is P_IMPLIED positive definite 
    # Expected value = TRUE
    
    improper_sol_ml[s, 2, n] = any(eigen(R_LVM_ML)$values<0)
    
    
    # Improper solution 3 : is variance of the residual positive
    # Expected value = TRUE
    
    improper_sol_ml[s, 3, n] = any(RESID_VAR_ML<0)
    
    # Improper solution 4 : number of out of range cor(y_jh, eta_j) 
    # Expected value : -1 <= cor(y_jh, eta_j) <=1
    
    improper_sol_ml[s, 4, n] = any(abs(STD_LAMBDA_ML)>1)
    
    # Improper solution 5 : is SIGMA_IMPLIED positive definite 
    # Expected value = TRUE
    
    improper_sol_ml[s, 5, n] = any(eigen(SIGMA_ML)$values<0)
    
    # Improper solution 6 : number of out of range R^2_i
    # Expected value : 0 <= R^2_i <=1
    
    improper_sol_ml[s, 6, n] = any(R2_ML< 0 | R2_ML>1)
    
    # Improper solution 7 : is PSI positive definite 
    # Expected value = TRUE
    
    improper_sol_ml[s, 7, n] = any(eigen(PSI_ML)$values<0) 
    }, silent = FALSE)

  }
}

dimnames(improper_sol_ml) = list(1:n_simu, c("RHO_JH", 
                                             "P_IMPLIED", 
                                             "THETA_JH", 
                                             "STD_LAMBDA", 
                                             "SIGMA_IMPLIED",
                                             "R2",
                                             "PSI"), N)

Table_improper = sapply(1:length(N), 
                        function(x) colMeans(improper_sol[, , x]))
colnames(Table_improper) = N
Table_improper[-c(1, 3), ]


Table_improper_ml = sapply(1:length(N), 
                        function(x) colMeans(improper_sol_ml[, , x]))
colnames(Table_improper_ml) = N
Table_improper_ml


#####################################
############ Table SIGMA ############
#####################################

Table1 = matrix(rowMeans(sigma_hat), 4, 3)
rownames(Table1) = c("Sigma_11", "Sigma_22", "Sigma_33", "Sigma_44")
colnames(Table1) = c("d_LS(empirical_S, Sigma)",
                     "d_LS(implied_S, Sigma)",
                     "d_LS(implied_S, empirical_S)")

Table1

#####################################
############ Table LAMBDA ###########
#####################################

Table2 = matrix(0, 18, 8)
Table2[, 1] = unlist(lambda)
Table2[, 2] = rowMeans(lambda_hat_svd)
Table2[, 3] = rowMeans(abs(lambda_hat_svd - Reduce("c", lambda)))
Table2[, 4] = rowMeans(lambda_hat_ml)
Table2[, 5] = rowMeans(abs(lambda_hat_ml - Reduce("c", lambda)))
Table2[, 6] = apply(lambda_hat_svd, 1, sd)
Table2[, 7] = apply(lambda_hat_ml, 1, sd)
Table2[, 8] = rowMeans(std_err_ml[1:18, ])

rownames(Table2) = paste("lambda_", rep(1:6, each = 3), rep(1:3,6), sep = "")
colnames(Table2) = c("Population", "Mean_SVD", "MAD_SVD",
                     "Mean_ML", "MAD_ML", "std_SVD", 
                     "std_ML", "std_RML")

round(Table2, 3)

#################################################
############ Table residual variance ############
#################################################

Table3 = matrix(0, 6, 8)
Table3[, 1] = 1-0.7^2
Table3[, 2] = rowMeans(var_hat_svd)
Table3[, 3] = rowMeans(abs(var_hat_svd - (1-0.7^2)))
Table3[, 4] = rowMeans(var_hat_ml)
Table3[, 5] = rowMeans(abs(var_hat_ml - (1-0.7^2)))
Table3[, 6] = apply(var_hat_svd, 1, sd)
Table3[, 7] = apply(var_hat_ml, 1, sd)
Table3[, 8] = rowMeans(std_err_ml[56:61, ])

rownames(Table3) = c("theta_51", "theta_52", "theta_53",
                     "theta_61", "theta_62", "theta_63")

colnames(Table3) = c("Population", "Mean_SVD", "MAD_SVD",
                     "Mean_ML", "MAD_ML", 
                     "std_SVD", "std_ML", "std_RML")
round(Table3, 3)

################################
#            Table 4           #
##########    OMEGA    #########
################################

Table4 = matrix(0, 12, 7)
Table4[, 1] = unlist(omega)
Table4[, 2] = rowMeans(omega_hat_svd)
Table4[, 3] = rowMeans(abs(omega_hat_svd - Reduce("c", omega)))
Table4[, 4] = rowMeans(omega_hat_ml)
Table4[, 5] = rowMeans(abs(omega_hat_ml - Reduce("c", omega)))
Table4[, 6] = apply(omega_hat_svd, 1, sd)
Table4[, 7] = apply(omega_hat_ml, 1, sd)

rownames(Table4) = paste("omega_", rep(1:4, each = 3), rep(1:3,4), sep = "")
colnames(Table4) = c("Population", "Mean_SVD", "MAD_SVD",  
                     "Mean_ML", "MAD_ML", "std_SVD", "std_ML")

round(Table4, 3)

#######################################################
####################### Table5 ########################
############ Table structural model parameters ########
#######################################################

Table5 = matrix(NA, 11, 8)
Table5[1:4, 1] = GAMMA[GAMMA!=0]
Table5[1:4, 2] = rowMeans(gamma_hat_svd)
Table5[1:4, 3] = rowMeans(abs(gamma_hat_svd - GAMMA[GAMMA!=0]))
Table5[1:4, 4] = rowMeans(gamma_hat_ml)
Table5[1:4, 5] = rowMeans(abs(gamma_hat_ml - GAMMA[GAMMA!=0]))
Table5[1:4, 6] = apply(gamma_hat_svd, 1, sd)
Table5[1:4, 7] = apply(gamma_hat_ml, 1, sd)
Table5[1:4, 8] = rowMeans(std_err_ml[25:28, ])

Table5[5:6, 1] = BETA[BETA!=0]
Table5[5:6, 2] = rowMeans(beta_hat_svd)
Table5[5:6, 3] = rowMeans(abs(beta_hat_svd - BETA[BETA!=0]))
Table5[5:6, 4] = rowMeans(beta_hat_ml)
Table5[5:6, 5] = rowMeans(abs(beta_hat_ml - BETA[BETA!=0]))
Table5[5:6, 6] = apply(beta_hat_svd, 1, sd)
Table5[5:6, 7] = apply(beta_hat_ml, 1, sd) 
Table5[5:6, 8] = rowMeans(std_err_ml[29:30, ])

Table5[7:9, 1] = PSI[as.vector(upper.tri(PSI, diag = TRUE))]
Table5[7:9, 2] = rowMeans(psi_hat_svd)
Table5[7:9, 3] = rowMeans(abs(psi_hat_svd - PSI[as.vector(upper.tri(PSI, diag = TRUE))]))
Table5[7:9, 4] = rowMeans(psi_hat_ml)
Table5[7:9, 5] = rowMeans(abs(psi_hat_ml - PSI[as.vector(upper.tri(PSI, diag = TRUE))]))
Table5[7:9, 6] = apply(psi_hat_svd, 1, sd)
Table5[7:9, 7] = apply(psi_hat_ml, 1, sd)


Table5[10:11, 1] = R2
Table5[10:11, 2] = rowMeans(r2_hat_svd)
Table5[10:11, 3] = rowMeans(abs(r2_hat_svd - R2))
Table5[10:11, 4] = rowMeans(r2_hat_ml)
Table5[10:11, 5] = rowMeans(abs(r2_hat_ml - R2))
Table5[10:11, 6] = apply(r2_hat_svd, 1, sd)
Table5[10:11, 7] = apply(r2_hat_ml, 1, sd)

rownames(Table5) = c("gamma_11", "gamma_12", "gamma_23", "gamma_24",
                     "beta_21", "beta_12", "psi11", "psi12", "psi22",
                     "R2_1", "R2_2")

colnames(Table5) = c("Population", "Mean_SVD", "MAD_SVD",
                     "Mean_ML", "MAD_ML", 
                     "std_SVD", "std_ML", "std_RML")

round(Table5, 3)

##################################
############# Table6 #############
############ Table rho ###########
##################################

Table6 = matrix(NA, 15, 8)
Table6[, 1] = R[as.vector(upper.tri(R))]
Table6[, 2] = rowMeans(rho_hat_svd)
Table6[, 3] = rowMeans(abs(rho_hat_svd - R[as.vector(upper.tri(R))]))
Table6[, 4] = rowMeans(rho_hat_ml)
Table6[, 5] = rowMeans(abs(rho_hat_ml - R[as.vector(upper.tri(R))]))
Table6[, 6] = apply(rho_hat_svd, 1, sd)
Table6[, 7] = apply(rho_hat_ml, 1, sd)
Table6[c(1:6, 15), 8] = rowMeans(std_err_ml[c(19:24, 31), ])


rownames(Table6) = c("p12", "p13", "p23", "p14", "p24",
                     "p34", "p15", "p25", "p35", "p45",
                     "p16", "p26", "p36", "p46", "p56")

colnames(Table6) = c("Population", "Mean_SVD", "MAD_SVD",
                     "Mean_ML", "MAD_ML", 
                     "std_err_SVD", "std_err_ML", "std_RML")

round(Table6, 3)

#####################################
############ Table 8 ################
# Empirical rejection probabilities #
#####################################

z = (gamma_hat_svd[2, ]-gamma_hat_svd[3, ])/sd(gamma_hat_svd[2, ]-gamma_hat_svd[3, ])

Table8 = matrix(NA, 1, 3)
Table8[1, 1] = mean(2*pnorm(abs(z), lower.tail = F) < 0.05)
Table8[1, 2] = mean(vcov_pval<0.05)
Table8[1, 3] = mean(lrt_pval<0.05)

colnames(Table8) = c("SVD_SEM", "ML-SEM_ZSCORE", "ML-SEM_LRT")
rownames(Table8) = c("% of false alarms")
round(Table8, 3)

#########################################
###      hypothesis testing with      ###
###           ML-SEM, eigenSEM        ###
#########################################

# Bootstrap testing
nsimu = 1000
nboot = 1000
N = c(300, 600, 1200)

decision_svd = matrix(NA, nsimu, length(N))
decision_ml = matrix(NA, nsimu, length(N))

T_LS_SVD = matrix(NA, nsimu, length(N))
Tb_LS_SVD = array(NA, c(nboot, nsimu, length(N)))
f_ml = matrix(NA, nsimu, length(N))

set.seed(20091979) #my date of birth

for(n in 1:length(N)){
  for (s in 1:nsimu){
    X = mvrnorm(N[n], rep(0, 18), SIGMA, empirical = FALSE)
    colnames(X) = paste("X", rep(1:6, each = 3), rep(1:3, 6), sep ="")
    S = cov(X)
    
    Y = list(X1 = X[, 1:3], X2 = X[, 4:6], X3 = X[, 7:9],
             X4 = X[, 10:12], X5 = X[, 13:15], X6 = X[, 16:18])
    
    C = matrix(c(0, 0, 0, 0, 1, 0,
                 0, 0, 0, 0, 1, 0,
                 0, 0, 0, 0, 0, 1,
                 0, 0, 0, 0, 0, 1,
                 0, 0, 0, 0, 0, 1,
                 0, 0, 0, 0, 1, 0), 6, 6, byrow = TRUE)
    
    colnames(C) = rownames(C) = names(Y)
    
    ##########
    #  SVD   #
    ##########  
    
    fit.svd = svdSEM(Y, C, scale = FALSE, 
                     mode = c(rep("formative", 4), 
                              rep("reflective", 2)), 
                     bias = FALSE)
    
    P_LVM_SVD = fit.svd$P_IMPLIED
    eigenval = eigen(P_LVM_SVD)$values
    
    if(!any(eigenval<0)){
      T_LS_SVD[s, n] = fit.svd$T_LS
      try(
        {
          
          ##########
          #   ML   #
          ##########
          
          init_ml_with_S =
            as.vector(c(Reduce("c", fit.svd$lambda),
                        P_LVM_SVD[1:4, 1:4][upper.tri( P_LVM_SVD[1:4, 1:4])],
                        fit.svd$gamma[1, 1:2], fit.svd$gamma[2, 3:4],
                        fit.svd$beta[1, 2], fit.svd$beta[2, 1],
                        P_LVM_SVD[5,6],
                        S[1:3, 1:3][upper.tri(S[1:3, 1:3], diag = TRUE)],
                        S[4:6, 4:6][upper.tri(S[4:6, 4:6], diag = TRUE)],
                        S[7:9, 7:9][upper.tri(S[7:9, 7:9], diag = TRUE)],
                        S[10:12, 10:12][upper.tri(S[10:12, 10:12], diag = TRUE)],
                        apply(X[, 13:15], 2, var)- fit.svd$lambda[[5]]^2,
                        apply(X[, 16:18], 2, var)- fit.svd$lambda[[6]]^2
            ))
          
          
          fit.ml = solnp(pars = init_ml_with_S, 
                         fun=F1, eqfun=heq1, 
                         eqB = rep(0,4), S = S, 
                         control = list(trace = 0))
          
          
          # ML output
          f_ml[s, n] = F1(fit.ml$pars, S)
        })
    }
    
    # transform the data sets in the way proposed
    # by Yuan & Hayashi (2003)
    
    if(!any(eigenval<0)){
      Z_SVD = scaleDataSet(X, fit.svd$SIGMA_IMPLIED)
      Z_SVD = lapply(split(data.frame(t(Z_SVD)),
                           as.factor(rep(1:6, each = 3))),
                     t)
    }

  
    ################################
    # Bootstrap hypotestis testing #
    ################################
    
    for (b in 1:nboot){
      ind = sample(1:N[n], replace = TRUE)
      
      ##########
      # SVD  #
      ##########
      
      if(!any(eigenval<0)){
        Zb_SVD = lapply(Z_SVD, function(x) x[ind, ])
        Sb_SVD = cov2(Reduce("cbind", Zb_SVD), bias = FALSE)
        
        fit.svd = svdSEM(Zb_SVD, C, scale = FALSE, 
                         mode = c(rep("formative", 4), rep("reflective", 2)), 
                         bias = FALSE)
        
        # SIGMA_SEM_SVD = fit.svd$SIGMA_IMPLIED
        Tb_LS_SVD[b, s, n] = fit.svd$T_LS
      }
      
    }
    decision_ml[s, n] = pchisq((N[n]-1)*f_ml[s, n], 114, lower.tail = FALSE)<=0.05
    decision_svd[s, n] = T_LS_SVD[s, n] > quantile(Tb_LS_SVD[, s, n], 0.95, na.rm = TRUE)
  }
}

######################
####### Table 7 ######
# global test-of-fit #
######################
######################

Table7 = matrix(NA, 2, 3)

Table7 = rbind(colMeans(decision_svd),
                    colMeans(decision_ml))

rownames(Table7) = c("SVD-SEM", "ML-SEM")
colnames(Table7) = c("N=300", "N=600", "N=1200")
round(Table7, 3)

############################
############################
##### CASE STUDY: ECSI #####
############################
############################

ECSI = as.data.frame(read_excel("SEMFC/data/mobil.xls"))/10
A = list(CUSTOMER_E = ECSI[, c("CUEX1", "CUEX2", "CUEX3")],
         PERC_QUAL  = ECSI[, c("PERQ1", "PERQ2", "PERQ3", "PERQ4", 
                               "PERQ5", "PERQ6", "PERQ7")],
         PERC_VALUE = ECSI[, c("PERV1", "PERV2")],
         CUSTOMER_S = ECSI[, c("CUSA1", "CUSA2", "CUSA3")],
         CUSTOMER_L = ECSI[, c("CUSL1", "CUSL2", "CUSL3")]
)

#######
# ML #
#######

A2 = data.frame(Reduce("cbind", A))

sem.model <-  '
# latent variable definitions
    eta1 =~ CUEX1+CUEX2+CUEX3
    eta2 =~ PERQ1+PERQ2+PERQ3+PERQ4+PERQ5+PERQ6+PERQ7
    eta3 =~ PERV1+PERV2
    eta4 =~ CUSA1+CUSA2+CUSA3
    eta5 =~ CUSL1+CUSL2+CUSL3

    # Regressions
    eta2 ~ eta1
    eta3 ~ eta1 + eta2
    eta4 ~ eta1 + eta2 + eta3
    eta5 ~ eta4'

fit.sem.ml <- sem(sem.model, data=A2, estimator = "ML")

###########
# SVD-SEM #
###########

C = matrix(c(0, 0, 0, 0, 0,
             1, 0, 0, 0, 0,
             1, 1, 0, 0, 0,
             1, 1, 1, 0, 0,
             0, 0, 0, 1, 0),
           5, 5, byrow = FALSE)

fit = svdSEM(A, C, scale = FALSE, 
             mode = rep("reflective", length(A)), 
             bias = FALSE)


T_LS_SVD = fit$T_LS
T_LS_SEM = d_LS(cov(A2), inspect(fit.sem.ml, "cov.ov"))


######################
#      Table10       #
# Boostrap inference #
######################
boot_out = svdSEM_infer(fit, B = 1000)
estimate = parameterEstimates(fit.sem.ml, standardized = TRUE)

Table10 = data.frame(boot_out$std_lambda,
                     std_loadings_sem = estimate[1:18, 10],
                     pval_std_loadings_sem = estimate[1:18, 7])

round(Table10, 3)

############
# Table 11 #
############

#The model-implied correlation matrix of the latent variables.
P_LVM_SEM = inspect(fit.sem.ml, what = "cor.lv")

Table11 = data.frame(
  rho_tilda = fit$Ptilde[lower.tri(fit$Ptilde, diag = FALSE)],
  rho_hat_svd = fit$P_IMPLIED[lower.tri(fit$P_IMPLIED, diag = FALSE)],
  rho_hat_sem = P_LVM_SEM[lower.tri(P_LVM_SEM, diag = FALSE)])

rownames(Table11) = c("p12", "p13", "p14", "p15",
                      "p23", "p24", "p25",
                      "p34", "p35",
                      "p45")

round(Table11, 3)

############
# Table12a #
############
# Estimate of structural equations and p-value
b_svd = boot_out$beta[, 1]
g_svd = boot_out$gamma[, 1]
est_svd = c(g_svd[1:2], b_svd[1], g_svd[3], b_svd[-1])
pval_svd = c(boot_out$gamma[1:2, 4], boot_out$beta[1, 4], 
             boot_out$gamma[3, 4], boot_out$beta[-1, 4])

Table12a = data.frame(estimate_svd = est_svd, 
                      pval = pval_svd, 
                      estimate_sem = estimate[19:25, 4],
                      pval = estimate[19:25, 7]) 

round(Table12a, 3)

############
# Table12b #
############
# R-squared 
Table12b = data.frame(
  R2_SVD = fit$R2, 
  R2_SEM = inspect(fit.sem.ml, "R-square")[c("eta2",  "eta3", "eta4", "eta5")])

rownames(Table12b) = c("Perceived quality",
                       "Perceived value",
                       "Customer satisfaction",
                       "Customer loyalty")

round(Table12b, 3)


########################
# Table14: Test-of-fit #
########################

# Bootstrap testing: global fit
# Tranforms the data sets in the way proposed 
# by Yuan & Hayashi (2003)

gof = svdSEM_gof(fit, B = 2000) #~half of improper solutions
hist(gof$Tb_LS) ; abline(v = gof$T_LS, col = "red")
pval_svd = mean(gof$T_LS <= gof$Tb_LS, na.rm = T)
pval_svd


############
# Table 14 #
############
Test_of_fit_SVD = data.frame(stat = gof$T_LS, pval = pval_svd)
Test_of_fit_SEM = data.frame(
  CHISQ = fitMeasures(fit.sem.ml)["chisq"],
  DF = fitMeasures(fit.sem.ml)["df"],
  PVAL = fitMeasures(fit.sem.ml)["pvalue"],
  T_LS_SEM,
  CHISQ_ON_DF = fitMeasures(fit.sem.ml)["chisq"]/fitMeasures(fit.sem.ml)["df"],
  RMSEA = fitMeasures(fit.sem.ml)["rmsea"],
  RMSEA.CI.LOWER = fitMeasures(fit.sem.ml)["rmsea.ci.lower"],
  RMSEA.CI.UPPER = fitMeasures(fit.sem.ml)["rmsea.ci.upper"], 
  CFI = fitMeasures(fit.sem.ml)["cfi"])

rownames(Test_of_fit_SVD) = "SVD-SEM"
Test_of_fit_SVD
rownames(Test_of_fit_SEM) = "ML-SEM"
Test_of_fit_SEM

