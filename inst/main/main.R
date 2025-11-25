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