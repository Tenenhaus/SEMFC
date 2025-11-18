



############################
############################
##### CASE STUDY: ECSI #####
############################
############################

source('R/SEMFC/sem_f_c.R')
source('data/data_ecsi.R')
library(lavaan)


#######
# ML #
#######

fit.sem.ml <- sem(sem.model.ecsi, data=A2, estimator = "ML")

###########
# SVD-SEM #
###########

model_ecsi <- SemFC$new(data=A, relation_matrix = C_ecsi, mode=mode_ecsi, scale=F, bias=F)
model_ecsi$fit_svd()


T_LS_SVD <- model_ecsi$parameters$T_LS
T_LS_SEM <- d_LS(cov(A2), inspect(fit.sem.ml, "cov.ov"))


######################
#      Table10       #
# Boostrap inference #
######################


model_ecsi$svd_infer(B = 1000)
estimate <- parameterEstimates(fit.sem.ml, standardized = TRUE)

Table10 <- data.frame(model_ecsi$infer_estimate$std_lambda,
                      std_loadings_sem = estimate[1:18, 11],
                      pval_std_loadings_sem = estimate[1:18, 7])


round(Table10, 3)



############
# Table 11 #
############

#The model-implied correlation matrix of the latent variables.
P_LVM_SEM <- inspect(fit.sem.ml, what = "cor.lv")

Table11 <- data.frame(
  rho_tilda = model_ecsi$parameters$Ptilde[lower.tri(model_ecsi$parameters$Ptilde, diag = FALSE)],
  rho_hat_svd = model_ecsi$parameters$P_IMPLIED[lower.tri(model_ecsi$parameters$P_IMPLIED, diag = FALSE)],
  rho_hat_sem = P_LVM_SEM[lower.tri(P_LVM_SEM, diag = FALSE)])

rownames(Table11) <- c("p12", "p13", "p14", "p15",
                       "p23", "p24", "p25",
                       "p34", "p35",
                       "p45")

round(Table11, 3)


############
# Table12a #
############
# Estimate of structural equations and p-value
b_svd <- model_ecsi$infer_estimate$beta[, 1]
g_svd <- model_ecsi$infer_estimate$gamma[, 1]
est_svd <- c(g_svd[1:2], b_svd[1], g_svd[3], b_svd[-1])
pval_svd <- c(model_ecsi$infer_estimate$gamma[1:2, 4], model_ecsi$infer_estimate$beta[1, 4],
              model_ecsi$infer_estimate$gamma[3, 4], model_ecsi$infer_estimate$beta[-1, 4])

Table12a <- data.frame(estimate_svd = est_svd,
                      pval = pval_svd,
                      estimate_sem = estimate[19:25, 11],
                      pval = estimate[19:25, 7])

round(Table12a, 3)

############
# Table12b #
############
# R-squared
Table12b <- data.frame(
  R2_SVD = model_ecsi$parameters$R2,
  R2_SEM = inspect(fit.sem.ml, "R-square")[c("PERC_QUAL",  "PERC_VALUE", "CUSTOMER_S", "CUSTOMER_L")])

rownames(Table12b) <- c("Perceived quality",
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


model_ecsi$estimator <- 'svd'
model_ecsi$get_gof(B = 2000)


hist(model_ecsi$gof$bollen_stine$Tb_LS) ; abline(v = model_ecsi$gof$bollen_stine$T_LS, col = "red")
pval_svd <- model_ecsi$gof$bollen_stine$pval
pval_svd


############
# Table 14 #
############
Test_of_fit_SVD <- data.frame(stat = model_ecsi$gof$bollen_stine$T_LS, pval = pval_svd)
Test_of_fit_SEM <- data.frame(
  CHISQ = fitMeasures(fit.sem.ml)["chisq"],
  DF = fitMeasures(fit.sem.ml)["df"],
  PVAL = fitMeasures(fit.sem.ml)["pvalue"],
  T_LS_SEM,
  CHISQ_ON_DF = fitMeasures(fit.sem.ml)["chisq"]/fitMeasures(fit.sem.ml)["df"],
  RMSEA = fitMeasures(fit.sem.ml)["rmsea"],
  RMSEA.CI.LOWER = fitMeasures(fit.sem.ml)["rmsea.ci.lower"],
  RMSEA.CI.UPPER = fitMeasures(fit.sem.ml)["rmsea.ci.upper"],
  CFI = fitMeasures(fit.sem.ml)["cfi"])

rownames(Test_of_fit_SVD) <- "SVD-SEM"
Test_of_fit_SVD
rownames(Test_of_fit_SEM) <- "ML-SEM"
Test_of_fit_SEM



