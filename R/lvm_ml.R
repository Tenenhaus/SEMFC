#########################################
#       return \Sigma(\theta) for       #
#         hessian computation           #
#########################################
# source("R/ml_sem/lvm_ml/get_loadings.R")
# source("R/ml_sem/lvm_ml/get_correlation_coeff.R")
# source("R/ml_sem/lvm_ml/get_path_coeff.R")
# source("R/ml_sem/lvm_ml/get_bdiag.R")
# source("R/utils/get_lengths_theta.R")


#' Compute Implied Covariance Matrix for Latent Variable Model
#'
#' This function computes the implied covariance matrix (SIGMA) for a structural
#' equation model with latent variables. It extracts model parameters from the
#' optimization vector `x` and returns either the upper triangular values of the
#' implied covariance matrix or the full model structure.
#'
#' @param x Numeric vector containing all model parameters (loadings, correlations,
#'   path coefficients, and variance/covariance parameters).
#' @param block_sizes Integer vector specifying the number of indicators in each block.
#' @param mode Character vector indicating the measurement mode for each block
#'   ("formative" or "reflective").
#' @param lengths_parameter Integer vector specifying the length of each parameter
#'   group in `x` (loadings, exogenous correlations, gamma, beta, endogenous correlations,
#'   variance/covariance).
#' @param which_exo_endo List containing indices and structure information for
#'   exogenous and endogenous latent variables (output from `ind_exo_endo()`).
#' @param jac Logical value. If TRUE, returns only the upper triangular values of
#'   the implied covariance matrix (for Jacobian computation). If FALSE, returns
#'   the complete model structure (default: TRUE).
#' @param varnames Optional list of character vectors containing variable names
#'   for each block (default: NULL).
#'
#' @return If `jac = TRUE`, returns a numeric vector of upper triangular values
#'   (including diagonal) of the implied covariance matrix.
#'   If `jac = FALSE`, returns a list containing:
#'   \item{lambda}{List of loading vectors for each block.}
#'   \item{beta}{Matrix of path coefficients between endogenous latent variables.}
#'   \item{gamma}{Matrix of path coefficients from exogenous to endogenous variables.}
#'   \item{psi}{Covariance matrix of structural disturbances.}
#'   \item{R2}{Vector of R-squared values for endogenous variables.}
#'   \item{residual_variance}{List of residual variances for reflective indicators.}
#'   \item{S_composites}{List of variance-covariance matrices for formative composites.}
#'   \item{omega}{List of composite weights (omega) for formative blocks.}
#'   \item{P_EXO}{Correlation matrix of exogenous latent variables.}
#'   \item{P_ENDO}{Correlation matrix of endogenous latent variables.}
#'   \item{R_LVM}{Full correlation matrix of all latent variables.}
#'   \item{SIGMA_IMPLIED}{Implied covariance matrix of observed variables.}
#'
#' @details
#' The function follows these steps:
#' 1. Extracts loadings from parameter vector `x`
#' 2. Constructs correlation matrices for exogenous (P_EXO) and endogenous (P_ENDO) variables
#' 3. Builds path coefficient matrices gamma (G) and beta (B)
#' 4. Computes structural disturbance covariance (PSI) and R-squared values
#' 5. Constructs variance/covariance blocks (BDIAG) for formative and reflective indicators
#' 6. Computes the implied covariance matrix using the LISREL equation:
#'   SIGMA = L * R * L' + BDIAG
#'
#' @export
lvm_ml <- function(x, block_sizes, mode, lengths_parameter, which_exo_endo, jac = TRUE, varnames = NULL){

  n <- which_exo_endo$ind_exo
  m <- which_exo_endo$ind_endo



  start_indices_in_x <- cumsum(c(1, head(lengths_parameter, -1)))

  ######################################################
  ####### mapping of the loadings from x ###############
  ######################################################

  loadings <- get_loadings(x, block_sizes = block_sizes)
  # naming loadings
  if (!is.null(varnames)){
    loadings <- mapply(function(x, nms) {names(x) <- nms; x}, loadings, varnames, SIMPLIFY = FALSE)
  }


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
  rownames(G) <- names(m)

  ##################################################################
  ####### mapping of the path coefficients matrix B from x #########
  ##################################################################

  B <- get_path_coeff(x,
                      list_linked_exo_endo = which_exo_endo$Ji,
                      exo_or_endo_variable = m,
                      initial_start_index = start_indices_in_x[4])

  rownames(B) <- names(m)

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
  # name the value
  residual_variance <- unname(split(
    `names<-`(
      unlist(residual_variance),
      paste0(".", names(unlist(unname(loadings[mode == "reflective"]))))
    ),
    rep(seq_along(residual_variance), lengths(residual_variance))
  ))


  # Get the variance matrices for reflective blocks
  S_composites <- BDIAG[mode == 'formative']

  ########################################################################
  ############## Get omega ###############################################
  ########################################################################

  omega <- mapply(function(Sjj, lambda_j) solve(Sjj) %*% lambda_j,
                  S_composites, loadings[mode == "formative"], SIMPLIFY = FALSE)

  for (b in seq_len(length(omega))){
    rownames(omega[[b]]) <- names(loadings[mode == "formative"][[b]])
  }



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
