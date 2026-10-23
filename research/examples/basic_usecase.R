
# exemple

library(cSEM)
library(devtools)
load_all()

set.seed(27)
source('research/model/model_mixed.R')
Y <- Y_2
X <- X_2

set.seed(27)
model <- SemFC$new(data=Y, relation_matrix = C, mode=mode, estimator = "svd")

model$fit(infer = T, B = 100)
model$fit(infer = F)
model$summary(standardized = T, effect = T, all_measures = T)
est = model$parameterEstimates(standardized = T)
gof = model$get_gof('all')
print(model$get_estimate('psi'))

modelml <- SemFC$new(data=Y, relation_matrix = C, mode=mode)
modelml$fit(infer=T, tol=1e-08)
modelml$fit(infer=F, tol=1e-10)
modelml$summary(standardized = T, effect = T, all_measures = T)
est = modelml$parameterEstimates(standardized = T)
gofml = modelml$get_gof('all')
print(modelml$get_estimate('psi'))



source('inst/model/model_ecsi.R')
model_ecsi <- SemFC$new(data=A, relation_matrix = C_ecsi, mode=mode_ecsi, estimator = "ml")
model_ecsi$fit(infer=T)
model_ecsi$summary(standardized = T, effect = F, all_measures = F)
print(model_ecsi$get_estimate('psi'))


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




data(ECSI)
ECSI = ECSI/10
A = list(CUSTOMER_E = ECSI[, c("CUEX1", "CUEX2", "CUEX3")],
        PERC_QUAL  = ECSI[, c("PERQ1", "PERQ2", "PERQ3",
                              "PERQ4", "PERQ5", "PERQ6", "PERQ7")],
        PERC_VALUE = ECSI[, c("PERV1", "PERV2")],
        CUSTOMER_S = ECSI[, c("CUSA1", "CUSA2", "CUSA3")],
        CUSTOMER_L = ECSI[, c("CUSL1", "CUSL2", "CUSL3")])

C = matrix(c(0, 0, 0, 0, 0,
            1, 0, 0, 0, 0,
            1, 1, 0, 0, 0,
            1, 1, 1, 0, 0,
            0, 0, 0, 1, 0), 5, 5,
            byrow = FALSE)

colnames(C) = rownames(C) = names(A)

# Initialize SemFC object with data, model specification
# and svd-based estimation method.

sem_model <- SemFC$new(data = A,
                      relation_matrix = C,
                      mode = rep("reflective", 5),
                      scale = FALSE,
                      estimator = "svd")

# Fit the full model with bootstrap inference and goodness-of-fit.

sem_model$fit(infer = TRUE, B = 1000)




