
#importing data
source('sep/sepdata.R')
# importing model
source('sep/sepmodel.R')
#importing lavaan
library(lavaan)
#importing semfc
source('R/SEMFC/sem_f_c.R')

X <- as.matrix(sep_df)
X <- scale(X)

Y <- list(clinique = X[, col_clinique], biologique = X[, col_biologique])

#lavaan

fit.sep.lav <- sem(sep.model, data=X, estimator = "ML", likelihood="wishart")
estimate <- parameterEstimates(fit.sep.lav, standardized = TRUE)

summary(fit.sep.lav, fit.measures = T)

#semfc

model <- SemFC$new(data=Y, relation_matrix = C, mode=mode, scale=F, bias=F)
model$fit_svd()
model$fit_ml(initialisation_svd = TRUE)

#bootstrap inference
model$svd_infer(500)
model$boot_svd$gof$pval
model$boot_svd$lambda
model$boot_svd$gamma

#ml inference
model$ml_infer()


model$reliability('svd', 'Dillon')
model$reliability_value

model$summary('svd')
model$summary('ml')