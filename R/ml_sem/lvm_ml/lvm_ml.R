#########################################
#       return \Sigma(\theta) for       #
#         hessian computation           #
#########################################
source("R/ml_sem/lvm_ml/get_loadings.R")
source("R/ml_sem/lvm_ml/get_correlation_coeff.R")
source("R/ml_sem/lvm_ml/get_path_coeff.R")
source("R/ml_sem/lvm_ml/get_bdiag.R")
source("R/utils/get_lengths_theta.R")



lvm_ml <- function(x, block_sizes, mode, lengths_parameter, which_exo_endo, jac = TRUE){

  n <- which_exo_endo$ind_exo
  m <- which_exo_endo$ind_endo



  start_indices_in_x <- cumsum(c(1, head(lengths_parameter, -1)))

  ######################################################
  ####### mapping of the loadings from x ###############
  ######################################################

  loadings <- get_loadings(x, block_sizes = block_sizes)

  ##################################################################
  ####### mapping of the exogeneous correlation matrix from x ######
  ##################################################################

  P_EXO <- get_correlation_coeff(x,
                                 latent_variables = n,
                                 start_index = start_indices_in_x[2])

  ##################################################################
  ####### mapping of the path coefficients matrix G from x #########
  ##################################################################

  G <- get_path_coeff(x,
                      list_linked_exo_endo = which_exo_endo$Hi,
                      exo_or_endo_variable = n,
                      initial_start_index = start_indices_in_x[3])
  rownames(G) = names(m)

  ##################################################################
  ####### mapping of the path coefficients matrix B from x #########
  ##################################################################

  B <- get_path_coeff(x,
                      list_linked_exo_endo = which_exo_endo$Ji,
                      exo_or_endo_variable = m,
                      initial_start_index = start_indices_in_x[4])

  rownames(B) = names(m)

  ##################################################################
  ####### mapping of the endogeneous correlation matrix from x #####
  ##################################################################

  P_ENDO <- get_correlation_coeff(x,
                                  latent_variables = m,
                                  start_index = start_indices_in_x[5])

  ##################################################################
  ########################## Computation of PSI, R2 ################
  ##################################################################

  PSI <-  (diag(NROW(B)) - B)%*%P_ENDO%*%t((diag(NROW(B)) - B)) - G%*%P_EXO%*%t(G)
  R2 <- 1-diag(PSI)

  ########################################################################
  ########  Correlations between Latent/emergent variables   #############
  ########################################################################

  R <- rbind(cbind(P_EXO, P_EXO%*%t(G)%*%t(solve(diag(NROW(B)) - B))),
            cbind(solve(diag(NROW(B)) - B)%*%G%*%P_EXO, P_ENDO))

  ########################################################################
  ############## Compute the variance blocks #############################
  ########################################################################

  BDIAG <- get_bdiag_bis(x,
                         mode = mode,
                         block_sizes = block_sizes,
                         initial_start_index_cov = start_indices_in_x[6])

  # Get the residual variance for reflective blocks
  residual_variance <- lapply(BDIAG[mode == 'reflective'], diag)
  # Get the variance matrices for reflective blocks
  S_composites <- BDIAG[mode == 'formative']

  ########################################################################
  ############## Get omega ###############################################
  ########################################################################

  omega <- mapply(function(Sjj, lambda_j) solve(Sjj) %*% lambda_j,
                  S_composites, loadings[mode == "formative"], SIMPLIFY = FALSE)

  ########################################################################
  ####### Compute the implied covariance matrix implied by the model #####
  ########################################################################
  formative_loadings <- lapply(loadings[mode == "formative"], function(lambda) { lambda %*% t(lambda) })
  BDIAG[mode == 'formative'] <- mapply("-", S_composites, formative_loadings, SIMPLIFY = FALSE)

  L <- as.matrix(Matrix::bdiag(loadings))
  BDIAG <- as.matrix(Matrix::bdiag(BDIAG))

  implied_S <- L%*%R%*%t(L) + BDIAG


  out <- list(
    lambda = loadings,
    beta = B,
    gamma = G,
    psi = PSI,
    R2 = R2,
    residual_variance = residual_variance,
    S_composites = S_composites,
    omega = omega,
    P_EXO = P_EXO,
    P_ENDO = P_ENDO,
    R_LVM = R,
    SIGMA_IMPLIED  = implied_S
  )

  if(jac){
    return(implied_S[upper.tri(implied_S, diag = T)])
  }else{
    return(out)
  }

}
