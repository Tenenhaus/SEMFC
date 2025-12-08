##################################################################
# Reproducible results for Tenenhaus, Tenenhaus & Dijkstra paper #
#              submitted to ADAC in July 15th, 2024              #
##################################################################

#############################################
# Remove all objects from the R environment #
#############################################
rm(list = ls())


source('inst/main/monte_carlo_sim.R')
source('inst/main/table_parameters.R')
source('inst/main/improper_solutions_simulation.R')
source('inst/main/hypothesis_testing.R')
source('inst/main/case_study_ecsi.R')




# exemple

library(cSEM)

source('inst/model/model_mixed.R')
Y <- Y_2
X <- X_2

set.seed(20091979)
model <- SemFC$new(data=Y, relation_matrix = C, mode=mode, scale=F, bias=F)
model$fit('svd', B = 2000)
model$summary()

modelml <- SemFC$new(data=Y, relation_matrix = C, mode=mode, scale=F, bias=F)
modelml$fit('ml')
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


