
# exemple

library(cSEM)
library(devtools)
load_all()

set.seed(27)
source('inst/model/model_mixed.R')
Y <- Y_2
X <- X_2

set.seed(27)
model <- SemFC$new(data=Y, relation_matrix = C, mode=mode, estimator = "svd")
# model$fit_svd()
# model$get_gof()
model$fit(infer = T, B = 100)
model$fit(infer = F)
model$summary(standardized = T, effect = T, all_measures = T)

modelml <- SemFC$new(data=Y, relation_matrix = C, mode=mode)
modelml$fit(infer=T, tol=1e-08)
modelml$summary(standardized = T, effect = T, all_measures = T)




source('inst/model/model_ecsi.R')
model_ecsi <- SemFC$new(data=A, relation_matrix = C_ecsi, mode=mode_ecsi, estimator = "ml")
model_ecsi$fit(infer=T)
model_ecsi$summary()



fit.csem <- csem(.data = X,
           .model = sem.model,
           .approach_weights = "PLS-PM",
           .PLS_weight_scheme_inner = "factorial",
           .approach_paths = "2SLS",
           .instruments = list( eta5 = c("eta1", "eta2", "eta3", "eta4"),
                                eta6 = c("eta1", "eta2", "eta3", "eta4")),
           .PLS_ignore_structural_model = TRUE, .tolerance = 1e-8,
           .disattenuate = TRUE)

summarize(fit.csem)