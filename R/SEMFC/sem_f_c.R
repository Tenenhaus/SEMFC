library(R6)



#import utils functions
source('R/utils/get_parameter_model_sem.R')
source('R/utils/ind_exo_endo.R')
source('R/utils/get_lengths_theta.R')
#import functions from svd module
source('R/svd_sem/svdSEM.R')
source('R/svd_sem/parameters_svd.R')
source('R/svd_sem/svdSEM_infer.R')
#import functions from ml module
source('R/ml_sem/F1.R')
source('R/ml_sem/mlSEM.R')
source('R/ml_sem/mlSEM_infer.R')

library(Matrix)
library(knitr)
library(pheatmap)

library(MASS)
library(lavaan)





SemFC <- R6Class(
  "SemFC",
  public = list(
    # Attributs
    data = NULL,
    relation_matrix = NULL,
    which_exo_endo = NULL,
    scale = FALSE,
    mode = NULL,
    cov_S = NULL,
    bias = FALSE,
    svd_result = NULL,
    n_blocks = NULL,
    n_row = NULL,
    block_sizes = NULL,
    lengths_theta = NULL,
    S_composites = NULL,
    boot_svd = NULL,
    ml_infer_estimate = NULL,
    SD = NULL,
    VCOV = NULL,
    svd_parameters = list(),
    ml_parameters = list(),







    # Méthode d'initialisation
    initialize = function(data, relation_matrix, mode,scale, bias) {
      self$data <- data
      self$relation_matrix <- relation_matrix
      self$scale <- ifelse(is.null(scale), FALSE, scale)
      self$bias <- ifelse(is.null(bias), FALSE, bias)
      self$mode <- mode

      parameter_model <- get_parameter_model_sem(data, mode)
      self$n_blocks <- parameter_model$n_blocks
      self$n_row <- parameter_model$n_row
      self$block_sizes <- parameter_model$block_sizes
      self$cov_S <- parameter_model$S
      self$S_composites <- parameter_model$S_diag_composites

      which_exo_endo <- ind_exo_endo(relation_matrix)
      self$which_exo_endo <-which_exo_endo

      self$lengths_theta <- get_lengths_theta(relation_matrix, self$block_sizes, mode)




    },



    # Méthode fit utilisant la technique SVD
    fit_svd = function() {
      svd_result <- svdSEM(self$data,
                           self$relation_matrix,
                           self$scale,
                           self$mode,
                           self$bias)

      self$svd_parameters <- svd_result
      theta_svd <- parameters_svd(lambda = svd_result$lambda,
                                  P_EXO = svd_result$P_EXO,
                                  G = svd_result$gamma,
                                  B = svd_result$beta,
                                  P_ENDO = svd_result$P_ENDO,
                                  residual_variance = svd_result$residual_variance,
                                  S_composites = self$S_composites,
                                  mode = self$mode)

      self$svd_parameters$theta <- theta_svd
      self$svd_parameters$F <- F1(theta_svd, self$cov_S, self$block_sizes, self$mode, self$lengths_theta, self$which_exo_endo)


    },
    svd_infer = function(B = 1000, verbose = TRUE){
      if (is.null(self$svd_parameters)) {
        self$fit_svd()
      }
      boot_out <- svdSEM_infer(self$svd_parameters, B = 1000, verbose = TRUE)
      self$boot_svd <- boot_out
    },


    fit_ml = function(initialisation_svd = TRUE) {

      block_sizes <- self$block_sizes
      C <- self$relation_matrix
      mode <- self$mode

      # Initialisation par SVD si demandé
      if (initialisation_svd) {
        if (is.null(self$svd_parameters$theta)) {
          self$fit_svd() }

        initial_params <- self$svd_parameters$theta
      } else {
        initial_params <- NULL
      }

      ml_sol <- mlSEM(initial_params, block_sizes, mode, self$cov_S, self$lengths_theta, self$which_exo_endo)
      theta_ml <- ml_sol$pars
      self$ml_parameters <- lvm_ml(x = theta_ml, block_sizes = block_sizes, mode =mode,
                                   lengths_parameter = self$lengths_theta, which_exo_endo = self$which_exo_endo, jac = F)
      self$ml_parameters$theta <- theta_ml
      self$ml_parameters$F <- F1(theta_ml, self$cov_S, self$block_sizes, self$mode, self$lengths_theta, self$which_exo_endo)


    },



    ml_infer = function(){
      theta_ml <- self$ml_parameters$theta
      block_sizes <- self$block_sizes
      mode <- self$mode
      S <- self$cov_S
      N <- self$n_row



      ml_infer_estimate <- mlSEM_infer(theta_ml, S, block_sizes, mode, self$lengths_theta, N, self$ml_parameters,self$which_exo_endo)

      self$ml_infer_estimate <- ml_infer_estimate
    },



    summary = function(){
      cat("\nCoefficients:\n")
      cat("lambda:\n")
      printCoefmat(model$ml_infer_estimate$estimate$lambda, P.values = TRUE, has.Pvalue = TRUE)
      cat("beta:\n")
      printCoefmat(model$ml_infer_estimate$estimate$beta, P.values = TRUE, has.Pvalue = TRUE)
      cat("gamma:\n")
      printCoefmat(model$ml_infer_estimate$estimate$gamma, P.values = TRUE, has.Pvalue = TRUE)
      cat("variance:\n")
      printCoefmat(model$ml_infer_estimate$estimate$residual_variance, P.values = TRUE, has.Pvalue = TRUE)


    }

  )
)



