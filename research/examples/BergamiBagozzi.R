library(devtools)
load_all()
library(lavaan)

source('inst/model/model_bergami.R')

sem_svd <- SemFC$new(data=A_bergami,
                     relation_matrix = C_bergami,
                     mode=mode_bergami,
                     scale = F, bias = T,
                     estimator = "svd")

sem_svd$fit(infer = T, B = 100)
sem_svd$summary()
est = sem_svd$parameterEstimates(standardized = T)


#########
# mlSEM #
#########

sem_ml <- SemFC$new(data=A_bergami,
                    relation_matrix = C_bergami,
                    mode=mode_bergami,
                    scale = F, bias = T,
                    estimator = "ml")

sem_ml$fit(infer = F)
sem_ml$summary(standardized = T, all_measures = T)
est = sem_ml$parameterEstimates(standardized = T)


##########
# lavaan #
##########



lavaan_ml <- sem(sem.model.bergami,
                 data=BergamiBagozzi2000,
                 estimator = "ML",
                 likelihood="wishart")


summary(lavaan_ml, standardized=TRUE, fit.measures=TRUE, rsquare=TRUE)
estimate = parameterEstimates(lavaan_ml, standardized=TRUE)



matrix(inspect(lavaan_ml, what = "cor.lv"), 3, 3)
sem_ml$estimate$P_IMPLIED

matrix(inspect(lavaan_ml, what = "cov.ov"), 11, 11)
sem_ml$estimate$SIGMA_IMPLIED




###############
# comparisons #
###############
inspect(lavaan_ml, what = "cor.lv")
round(sem_svd$estimate$P_IMPLIED, 3)
round(sem_ml$estimate$P_IMPLIED, 3)
matrix(inspect(lavaan_ml, what = "cor.lv"), nrow=5, ncol=5)


