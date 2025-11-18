


#####################################
############ Table SIGMA ############
#####################################

Table1 <- matrix(rowMeans(sigma_hat), 4, 3)
rownames(Table1) <- c("Sigma_11", "Sigma_22", "Sigma_33", "Sigma_44")
colnames(Table1) <- c("d_LS(empirical_S, Sigma)",
                     "d_LS(implied_S, Sigma)",
                     "d_LS(implied_S, empirical_S)")

Table1




#####################################
############ Table LAMBDA ###########
#####################################

Table2 <- matrix(0, 18, 8)
Table2[, 1] <- unlist(lambda)
Table2[, 2] <- rowMeans(lambda_hat_svd)
Table2[, 3] <- rowMeans(abs(lambda_hat_svd - Reduce("c", lambda)))
Table2[, 4] <- rowMeans(lambda_hat_ml)
Table2[, 5] <- rowMeans(abs(lambda_hat_ml - Reduce("c", lambda)))
Table2[, 6] <- apply(lambda_hat_svd, 1, sd)
Table2[, 7] <- apply(lambda_hat_ml, 1, sd)
Table2[, 8] <- rowMeans(std_err_ml[1:18, ])

rownames(Table2) <- paste0("lambda_", rep(1:6, each <- 3), rep(1:3, 6))
colnames(Table2) <- c("Population", "Mean_SVD", "MAD_SVD",
                     "Mean_ML", "MAD_ML", "std_SVD",
                     "std_ML", "std_RML")

round(Table2, 3)


#################################################
############ Table residual variance ############
#################################################

Table3 <- matrix(0, 6, 8)
Table3[, 1] <- 1-0.7^2
Table3[, 2] <- rowMeans(var_hat_svd)
Table3[, 3] <- rowMeans(abs(var_hat_svd - (1-0.7^2)))
Table3[, 4] <- rowMeans(var_hat_ml)
Table3[, 5] <- rowMeans(abs(var_hat_ml - (1-0.7^2)))
Table3[, 6] <- apply(var_hat_svd, 1, sd)
Table3[, 7] <- apply(var_hat_ml, 1, sd)
Table3[, 8] <- rowMeans(std_err_ml[56:61, ])

rownames(Table3) <- c("theta_51", "theta_52", "theta_53",
                     "theta_61", "theta_62", "theta_63")

colnames(Table3) <- c("Population", "Mean_SVD", "MAD_SVD",
                     "Mean_ML", "MAD_ML",
                     "std_SVD", "std_ML", "std_RML")
round(Table3, 3)



################################
#            Table 4           #
##########    OMEGA    #########
################################

Table4 <- matrix(0, 12, 7)
Table4[, 1] <- unlist(omega)
Table4[, 2] <- rowMeans(omega_hat_svd)
Table4[, 3] <- rowMeans(abs(omega_hat_svd - Reduce("c", omega)))
Table4[, 4] <- rowMeans(omega_hat_ml)
Table4[, 5] <- rowMeans(abs(omega_hat_ml - Reduce("c", omega)))
Table4[, 6] <- apply(omega_hat_svd, 1, sd)
Table4[, 7] <- apply(omega_hat_ml, 1, sd)

rownames(Table4) <- paste0("omega_", rep(1:4, each <- 3), rep(1:3, 4))
colnames(Table4) <- c("Population", "Mean_SVD", "MAD_SVD",
                     "Mean_ML", "MAD_ML", "std_SVD", "std_ML")

round(Table4, 3)



#######################################################
####################### Table5 ########################
############ Table structural model parameters ########
#######################################################

Table5 <- matrix(NA, 11, 8)
Table5[1:4, 1] <- GAMMA[GAMMA!=0]
Table5[1:4, 2] <- rowMeans(gamma_hat_svd)
Table5[1:4, 3] <- rowMeans(abs(gamma_hat_svd - GAMMA[GAMMA!=0]))
Table5[1:4, 4] <- rowMeans(gamma_hat_ml)
Table5[1:4, 5] <- rowMeans(abs(gamma_hat_ml - GAMMA[GAMMA!=0]))
Table5[1:4, 6] <- apply(gamma_hat_svd, 1, sd)
Table5[1:4, 7] <- apply(gamma_hat_ml, 1, sd)
Table5[1:4, 8] <- rowMeans(std_err_ml[25:28, ])

Table5[5:6, 1] <- BETA[BETA!=0]
Table5[5:6, 2] <- rowMeans(beta_hat_svd)
Table5[5:6, 3] <- rowMeans(abs(beta_hat_svd - BETA[BETA!=0]))
Table5[5:6, 4] <- rowMeans(beta_hat_ml)
Table5[5:6, 5] <- rowMeans(abs(beta_hat_ml - BETA[BETA!=0]))
Table5[5:6, 6] <- apply(beta_hat_svd, 1, sd)
Table5[5:6, 7] <- apply(beta_hat_ml, 1, sd)
Table5[5:6, 8] <- rowMeans(std_err_ml[29:30, ])

Table5[7:9, 1] <- PSI[as.vector(upper.tri(PSI, diag = TRUE))]
Table5[7:9, 2] <- rowMeans(psi_hat_svd)
Table5[7:9, 3] <- rowMeans(abs(psi_hat_svd - PSI[as.vector(upper.tri(PSI, diag = TRUE))]))
Table5[7:9, 4] <- rowMeans(psi_hat_ml)
Table5[7:9, 5] <- rowMeans(abs(psi_hat_ml - PSI[as.vector(upper.tri(PSI, diag = TRUE))]))
Table5[7:9, 6] <- apply(psi_hat_svd, 1, sd)
Table5[7:9, 7] <- apply(psi_hat_ml, 1, sd)


Table5[10:11, 1] <- R2
Table5[10:11, 2] <- rowMeans(r2_hat_svd)
Table5[10:11, 3] <- rowMeans(abs(r2_hat_svd - R2))
Table5[10:11, 4] <- rowMeans(r2_hat_ml)
Table5[10:11, 5] <- rowMeans(abs(r2_hat_ml - R2))
Table5[10:11, 6] <- apply(r2_hat_svd, 1, sd)
Table5[10:11, 7] <- apply(r2_hat_ml, 1, sd)

rownames(Table5) <- c("gamma_11", "gamma_12", "gamma_23", "gamma_24",
                     "beta_21", "beta_12", "psi11", "psi12", "psi22",
                     "R2_1", "R2_2")

colnames(Table5) <- c("Population", "Mean_SVD", "MAD_SVD",
                     "Mean_ML", "MAD_ML",
                     "std_SVD", "std_ML", "std_RML")

round(Table5, 3)


##################################
############# Table6 #############
############ Table rho ###########
##################################

Table6 <- matrix(NA, 15, 8)
Table6[, 1] <- R[as.vector(upper.tri(R))]
Table6[, 2] <- rowMeans(rho_hat_svd)
Table6[, 3] <- rowMeans(abs(rho_hat_svd - R[as.vector(upper.tri(R))]))
Table6[, 4] <- rowMeans(rho_hat_ml)
Table6[, 5] <- rowMeans(abs(rho_hat_ml - R[as.vector(upper.tri(R))]))
Table6[, 6] <- apply(rho_hat_svd, 1, sd)
Table6[, 7] <- apply(rho_hat_ml, 1, sd)
Table6[c(1:6, 15), 8] <- rowMeans(std_err_ml[c(19:24, 31), ])


rownames(Table6) <- c("p12", "p13", "p23", "p14", "p24",
                     "p34", "p15", "p25", "p35", "p45",
                     "p16", "p26", "p36", "p46", "p56")

colnames(Table6) <- c("Population", "Mean_SVD", "MAD_SVD",
                     "Mean_ML", "MAD_ML",
                     "std_err_SVD", "std_err_ML", "std_RML")

round(Table6, 3)

#####################################
############ Table 8 ################
# Empirical rejection probabilities #
#####################################

z <- (gamma_hat_svd[2, ]-gamma_hat_svd[3, ])/sd(gamma_hat_svd[2, ]-gamma_hat_svd[3, ])

Table8 <- matrix(NA, 1, 3)
Table8[1, 1] <- mean(2*pnorm(abs(z), lower.tail = F) < 0.05)
Table8[1, 2] <- mean(vcov_pval<0.05)
Table8[1, 3] <- mean(lrt_pval<0.05)

colnames(Table8) <- c("SVD_SEM", "ML-SEM_ZSCORE", "ML-SEM_LRT")
rownames(Table8) <- "% of false alarms"
round(Table8, 3)