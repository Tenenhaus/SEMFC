
library(devtools)
load_all()
library(lavaan)

source('inst/model/model_itflex.R')



sem_svd <- SemFC$new(data=A_itflex,
                     relation_matrix = C_itflex,
                     mode=mode_itflex,
                     scale = F, bias = T,
                     estimator = "svd")

sem_svd$fit(infer = T, B = 100)
sem_svd$summary(effect = T, all_measures = T)
est = sem_svd$parameterEstimates(standardized = T)

#########
# mlSEM #
#########

sem_ml <- SemFC$new(data = A_itflex,
                    relation_matrix = C_itflex,
                    mode = mode_itflex,
                    scale = F, bias = F,
                    estimator = "ml")

sem_ml$fit(infer = T)
sem_ml$summary()
est = sem_ml$parameterEstimates(standardized = T)
