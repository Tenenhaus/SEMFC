
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
model$fit(estimator = 'svd', B=2000)
model$summary()


modelml <- SemFC$new(data=Y, relation_matrix = C, mode=mode, scale=F, bias=F)
modelml$fit(estimator = 'ml', initialisation_svd  = TRUE)
modelml$summary()





