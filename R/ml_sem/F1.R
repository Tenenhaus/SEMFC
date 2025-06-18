########################################
# Objective function for ML estimation #
########################################
source("R/ml_sem/lvm_ml/lvm_ml.R")



F1 <- function(x, S, block_sizes, mode, lengths_parameter, which_exo_endo){


  implied_S <- lvm_ml(x, block_sizes, mode, lengths_parameter, which_exo_endo, jac = FALSE)$SIGMA_IMPLIED

  ########################################################################
  ###################### Compute log-likelihood  #########################
  ########################################################################

  opt <- log(det(implied_S)) +
    sum(diag(S%*%solve(implied_S))) -
    log(det(S)) -
    NCOL(S)

  return(opt)

}