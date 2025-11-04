library(R6)



#import utils functions
source('R/utils/get_parameter_model_sem.R')
source('R/utils/ind_exo_endo.R')
source('R/utils/get_lengths_theta.R')
source('R/utils/reliability.R')
source('R/utils/chi2sem.R')
#import functions from svd module
source('R/svd_sem/svdSEM.R')
source('R/svd_sem/parameters_svd.R')
source('R/svd_sem/svdSEM_infer.R')
source("R/svd_sem/svdSEM_gof.R")
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
    estimator = NULL,
    relation_matrix = NULL,
    which_exo_endo = NULL,
    scale = FALSE,
    mode = NULL,
    cov_S = NULL,
    bias = FALSE,
    svd_result = NULL,
    n_blocks = NULL,
    n_row = NULL,
    varnames = NULL,
    block_sizes = NULL,
    lengths_theta = NULL,
    S_composites = NULL,
    infer_estimate = NULL,
    boot_rep = NULL,
    reliability_value = NULL,
    SD = NULL,
    VCOV = NULL,
    gof = NULL,
    parameters = list(),







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
      self$varnames <- parameter_model$varnames
      self$block_sizes <- parameter_model$block_sizes
      self$cov_S <- parameter_model$S
      self$S_composites <- parameter_model$S_diag_composites

      which_exo_endo <- ind_exo_endo(relation_matrix)
      self$which_exo_endo <- which_exo_endo

      self$lengths_theta <- get_lengths_theta(self$which_exo_endo, self$block_sizes, self$mode)


    },



    # Méthode fit utilisant la technique SVD
    fit_svd = function() {
      svd_result <- svdSEM(self$data,
                           self$relation_matrix,
                           self$scale,
                           self$mode,
                           self$bias)

      self$parameters <- svd_result
      theta_svd <- parameters_svd(lambda = svd_result$lambda,
                                  P_EXO = svd_result$P_EXO,
                                  G = svd_result$gamma,
                                  B = svd_result$beta,
                                  P_ENDO = svd_result$P_ENDO,
                                  residual_variance = svd_result$residual_variance,
                                  S_composites = self$S_composites,
                                  mode = self$mode)

      self$parameters$theta <- theta_svd
      self$parameters$F <- F1(theta_svd, self$cov_S, self$block_sizes, self$mode, self$lengths_theta, self$which_exo_endo)
    },

    svd_infer = function(B = 1000, verbose = TRUE){
      if (is.null(self$parameters)) {
        self$fit_svd()
      }
      boot_out <- svdSEM_infer(self$parameters, B, verbose = TRUE)
      self$infer_estimate <- boot_out


    },


    fit_ml = function(initialisation_svd = TRUE) {

      block_sizes <- self$block_sizes
      mode <- self$mode

      # Initialisation par SVD si demandé
      if (initialisation_svd) {
        if (is.null(self$parameters$theta)) {
          self$fit_svd() }

        initial_params <- self$parameters$theta
      } else {
        len_theta <- sum(self$lengths_theta)
        initial_params <- runif(len_theta)
      }

      ml_sol <- mlSEM(initial_params, block_sizes, mode, self$cov_S, self$lengths_theta, self$which_exo_endo)
      theta_ml <- ml_sol$pars
      self$parameters <- lvm_ml(x = theta_ml, block_sizes = block_sizes, mode =mode,
                                   lengths_parameter = self$lengths_theta, which_exo_endo = self$which_exo_endo,
                                   jac = F, varnames = self$varnames)
      self$parameters$theta <- theta_ml
      self$parameters$F <- F1(theta_ml, self$cov_S, self$block_sizes, self$mode, self$lengths_theta, self$which_exo_endo)


    },



    ml_infer = function(){
      theta_ml <- self$parameters$theta
      block_sizes <- self$block_sizes
      mode <- self$mode
      S <- self$cov_S
      N <- self$n_row

      ml_infer_estimate <- mlSEM_infer(theta_ml, S, block_sizes, mode, self$lengths_theta, N, self$parameters,self$which_exo_endo)

      self$infer_estimate <- ml_infer_estimate$estimate
    },

    reliability = function(fit, metric='Dillon'){
      if (fit=='svd'){
        lambdas <- self$parameters$lambda
        residual_variances <- self$parameters$residual_variance
      }
      else if (fit=='ml'){
        lambdas <- self$parameters$lambda
        residual_variances <- self$parameters$residual_variance

      }

      res_reliability <- reliability(metric, lambdas, residual_variances)

      self$reliability_value <- res_reliability

    },

    get_gof = function(B = 1000){

      estimator <- self$estimator

      p <- sum(self$block_sizes)
      q <- sum(self$lengths_theta)
      F <- self$parameters$F
      N <- self$n_row
      S <- self$cov_S
      Sigma <- self$parameters$SIGMA_IMPLIED

      res_gof <- list()

      chi2 <- chi2sem(p, q, F, N)
      res_gof$chi2 <- chi2

      # basline test
      S_baseline <- diag(diag(S))
      F_baseline <- log(det(S_baseline)) + sum(diag(S%*%solve(S_baseline))) - log(det(S)) - NCOL(S)
      baseline <- chi2sem(p, p, F_baseline, N)
      res_gof$baseline <- baseline

      if (estimator == 'svd'){
        bollen_stine <- svdSEM_gof(self$parameters, B)
        res_gof$bollen_stine <- bollen_stine
      }


      cfi <- 1 - (max(chi2$test - chi2$df, 0)) /
           (max(baseline$test - baseline$df, chi2$test - chi2$df, 0))

      tli <- ( (baseline$test /  baseline$df) - (chi2$test / chi2$df) ) /
             ( (baseline$test /  baseline$df) - 1 )

      res_gof$cfi <- cfi
      res_gof$tli <- tli

      # loglik
      loglik_H0 <- -(N/2)*(p*log(2*pi) + log(det(Sigma)) + sum(diag(solve(Sigma) %*% S)))
      loglik_H1 <- -(N/2)*(p*log(2*pi) + log(det(S)) + sum(diag(solve(S) %*% S)))
      AIC <- -2 * loglik_H0 + 2 * q
      BIC <- -2 * loglik_H0 + q * log(N)
      SABIC <- -2 * loglik_H0 + q * log((N + 2) / 24)



      self$gof <- res_gof

    },



    fit = function(estimator, B = 1000, initialisation_svd = TRUE){
      self$estimator <- estimator
      self$boot_rep <- B
      if (estimator == 'svd'){
        self$fit_svd()
        self$svd_infer(B)
      } else if(estimator == 'ml'){
        self$fit_ml(initialisation_svd)
        self$ml_infer()
      }
      self$get_gof(B)

    },

    summary = function(){

      estimator <- self$estimator
      infer_estimate <- self$infer_estimate

      #gof
      # user test
      testchi2 <- self$gof$chi2$test
      dfchi2 <- self$gof$chi2$df
      pvalchi2 <- self$gof$chi2$pval

      # baseline test

      testbaseline <- self$gof$baseline$test
      dfbaseline <- self$gof$baseline$df
      pvalbaseline <- self$gof$baseline$pval

      #  vs
      cfi <- self$gof$cfi
      tli <- self$gof$tli




      # inference estimation
      lambda_infer <- infer_estimate$lambda
      beta_infer <- infer_estimate$beta
      gamma_infer <- infer_estimate$gamma
      residualvariance_infer <- infer_estimate$residual_variance


      cat("\n")
      cat("Estimator", sprintf("%44s\n", toupper(estimator)))
      cat("Number of model parameters", sprintf("%23d\n", sum(self$lengths_theta)))
      cat("\n")
      cat("Number of observations", sprintf("%32d\n", self$n_row))
      cat("\n\n")

      cat("Model Test User Model :\n")
      cat("                                                      \n")
      cat(sprintf("  Test statistic%38.3f\n", testchi2))
      cat(sprintf("  Degrees of freedom%33d\n", dfchi2))
      cat(sprintf("  P-value (Chi-square)%29.3f\n", pvalchi2))
      cat("\n")


      cat("Model Test Baseline Model :\n")
      cat("                                                      \n")
      cat(sprintf("  Test statistic%38.3f\n", testbaseline))
      cat(sprintf("  Degrees of freedom%33d\n", dfbaseline))
      cat(sprintf("  P-value %29.3f\n", pvalbaseline))
      cat("\n")


      cat("User Model versus Baseline Model:\n")
      cat("\n")
      cat(sprintf("  Comparative Fit Index (CFI)%22.3f\n", cfi))
      cat(sprintf("  Tucker-Lewis Index (TLI)%25.3f\n", tli))
      cat("\n")



      if (estimator == 'svd'){
        B <- self$boot_rep
        pvalbs <- self$gof$bollen_stine$pval

        cat("Bootstrap Test (Bollen Stine):\n")
        cat("\n")
        cat(sprintf("  Number of bootstrap replications%17d\n", B))
        cat(sprintf("  Bollen Stine bootstrap p-value%20.3f\n", pvalbs))
        cat("\n")
      }


      cat("\nCoefficients:\n")
      cat("lambda:\n")
      if (nrow(lambda_infer) != 0){
        printCoefmat(lambda_infer, P.values = TRUE, has.Pvalue = TRUE)
      }

      if (nrow(beta_infer) != 0){
        cat("beta:\n")
        printCoefmat(beta_infer, P.values = TRUE, has.Pvalue = TRUE)
      }
      cat("gamma:\n")
      if (nrow(gamma_infer) != 0){
        printCoefmat(gamma_infer, P.values = TRUE, has.Pvalue = TRUE)
      }

      cat("variance:\n")
      if (nrow(residualvariance_infer) != 0){
        printCoefmat(residualvariance_infer, P.values = TRUE, has.Pvalue = TRUE)
      }

    }

  )
)



