
# exemple

library(cSEM)

source('inst/model/model_mixed.R')
Y <- Y_2
X <- X_2

set.seed(1)
model <- SemFC$new(data=Y, relation_matrix = C, mode=mode, estimator = "svd")
model$fit_svd()
model$get_gof()
model$fit(infer = T, B = 100)
model$fit(infer = F)
model$summary()

modelml <- SemFC$new(data=Y, relation_matrix = C, mode=mode, scale=F, bias=F)
modelml$fit(infer=T)
modelml$summary()


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