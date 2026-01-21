
#' @import R6
#' @title SemFC Class
#'
#' @description
#' R6 class for estimating and analyzing Structural Equation Models (SEM) that
#' incorporate both latent factors and composite variables. Supports SVD-based
#' and maximum likelihood estimation methods.
#' @details
#' This class provides a complete framework for SEM analysis including:
#' \itemize{
#'   \item Model estimation using SVD or ML
#'   \item Statistical inference via bootstrap or asymptotic methods
#'   \item Goodness-of-fit assessment
#'   \item Reliability coefficients calculation
#'   \item Support for formative and reflective measurement models
#' }
#' @field estimator Character string specifying estimation method ("svd" or "ml")
#' @field data List containing data-related components:
#'   \itemize{
#'     \item \code{data}: The input data (list of blocks)
#'     \item \code{n_row}: Number of observations
#'     \item \code{cov_S}: Covariance matrix of observed variables
#'     \item \code{S_diag_composites}: List of covariance matrices for formative blocks
#'   }
#' @field model List containing model specification and parameters:
#'   \itemize{
#'     \item \code{relation_matrix}: Square matrix defining structural relationships
#'     \item \code{mode}: Character vector specifying measurement model types
#'     \item \code{n_blocks}: Number of measurement blocks
#'     \item \code{varnames}: List of variable names for each block
#'     \item \code{block_sizes}: Vector of sizes for each measurement block
#'     \item \code{dag}: Logical indicating if structural model is recursive
#'     \item \code{which_exo_endo}: List identifying exogenous/endogenous variables
#'     \item \code{lengths_theta}: Vector of parameter counts
#'     \item \code{p}: Total number of observed variables
#'     \item \code{q}: Total number of free parameters
#'     \item \code{r}: Number of formative blocks
#'     \item \code{dof}: Degrees of freedom
#'     \item \code{scale}: Logical indicating whether to standardize data
#'     \item \code{bias}: Logical indicating bias correction in covariance estimation
#'   }
#' @field estimate List containing all estimated model parameters
#' @field infer_estimate Data frame containing inference results (estimates, SE, z-values, p-values)
#' @field boot_rep Integer number of bootstrap replications
#' @field gof List containing goodness-of-fit statistics
#'
#' @examples
#' \dontrun{
#' # Create relation matrix
#' rel_matrix <- matrix(0, 3, 3)
#' rel_matrix[1, 3] <- 1
#' rel_matrix[2, 3] <- 1
#'
#' # Fit model
#' model <- SemFC$new(
#'   data = my_data,
#'   relation_matrix = rel_matrix,
#'   mode = c("reflective", "reflective", "reflective"),
#'   scale = TRUE,
#'   bias = FALSE
#' )
#' model$fit(estimator = "svd", B = 1000)
#' model$summary()
#' }
#' @export
SemFC <- R6Class(
  "SemFC",
  public = list(
    # Attributs
    estimator = NULL,
    data = list(),
    model = list(),
    estimate = list(),
    infer_estimate = NULL,
    boot_rep = NULL,
    gof = NULL,



    #' @description
    #' Create a new SemFC object and initialize model parameters
    #'
    #' @param data Data frame or matrix where each column represents an observed variable
    #' @param relation_matrix Square adjacency matrix (n_blocks x n_blocks) defining structural
    #'   paths between latent variables (1 = path exists, 0 = no path)
    #' @param mode Character vector of length n_blocks specifying measurement model type
    #'   ("formative" or "reflective") for each block
    #' @param estimator Character string specifying estimation method: "svd" or "ml"
    #' @param scale Logical indicating whether to standardize input data (default: FALSE)
    #' @param bias Logical indicating whether to apply bias correction in covariance
    #'   estimation (default: FALSE)
    #'
    #' @return A new `SemFC` object

    initialize = function(data, relation_matrix, mode, estimator = 'ml', scale = FALSE, bias = FALSE) {

      self$estimator <- estimator
      init <- get_parameter_model_sem(data, mode, relation_matrix, bias)
      self$model <- init$model
      self$data <- init$data

      self$model$scale <- scale
      self$model$bias <- bias

    },


    #' @description
    #' Fit the model using Singular Value Decomposition (SVD) method
    #'
    #' @details
    #' Estimates model parameters using SVD-based approach which is computationally
    #' efficient and provides good starting values for ML estimation.
    #'
    #' @return Invisible self (for method chaining)
    # Méthode fit utilisant la technique SVD
    fit_svd = function() {
      self$estimator <- 'svd'
      svd_result <- svdSEM(self$data$data,
                           self$model$relation_matrix,
                           self$model$scale,
                           self$model$mode,
                           self$model$bias)

      self$estimate <- svd_result
      theta_svd <- parameters_svd(lambda = svd_result$lambda,
                                  P_EXO = svd_result$P_EXO,
                                  G = svd_result$gamma,
                                  B = svd_result$beta,
                                  P_ENDO = svd_result$P_ENDO,
                                  residual_variance = svd_result$residual_variance,
                                  S_composites = self$data$S_diag_composites,
                                  model = self$model)

      self$estimate$theta <- theta_svd
      self$estimate$effect <- compute_effect(self$estimate$beta, self$estimate$gamma)
      self$gof$F <- F1(theta_svd, self$data$cov_S, self$model)
    },


     #' @description
    #' Perform statistical inference for SVD estimates using bootstrap
    #'
    #' @param B Integer number of bootstrap replications (default: 1000)
    #' @param verbose Logical indicating whether to print progress messages (default: TRUE)
    #'
    #' @details
    #' Uses non-parametric bootstrap to estimate standard errors, confidence intervals,
    #' and p-values for all model parameters.
    #'
    #' @return Invisible self (for method chaining)
    svd_infer = function(B = 1000, verbose = TRUE){
      if (is.null(self$estimate)) {
        self$fit_svd()
      }
      boot_out <- svdSEM_infer(self$estimate, B, verbose = TRUE)
      self$infer_estimate <- boot_out


    },

    #' @description
    #' Fit the model using Maximum Likelihood (ML) estimation
    #'
    #' @param initialisation_svd Logical indicating whether to use SVD estimates as
    #'   starting values (default: TRUE). If FALSE, random starting values are used.
    #' @param tol Numeric tolerance for convergence in optimization (default: 1e-8)
    #'
    #' @details
    #' Uses numerical optimization (via SOLNP) to minimize the ML fit function.
    #' SVD initialization is recommended for better convergence.
    #'
    #' @return Invisible self (for method chaining)
    fit_ml = function(initialisation_svd = TRUE, tol) {


      # Initialisation par SVD si demandé
      if (initialisation_svd) {
        self$fit_svd()
        self$estimator <- 'ml'
        initial_params <- self$estimate$theta
      } else {
        len_theta <- sum(self$model$lengths_theta)
        initial_params <- runif(len_theta)
      }

      ml_sol <- mlSEM(initial_params, self$data$cov_S, self$model, tol)
      theta_ml <- ml_sol$pars
      self$estimate <- lvm_ml(x = theta_ml, model = self$model, jac = F)
      self$estimate$T_LS <- d_LS(self$data$cov_S, self$estimate$SIGMA_IMPLIED)

      var_MVs <- lapply(self$data$data, function(x) diag(cov2(x, bias = self$model$bias)))
      std_lambda <- mapply("/", self$estimate$lambda, lapply(var_MVs, sqrt),  SIMPLIFY = FALSE)
      self$estimate$std_lambda <- std_lambda

      self$estimate$theta <- theta_ml
      self$estimate$effect <- compute_effect(self$estimate$beta, self$estimate$gamma)
      self$gof$F <- F1(theta_ml, self$data$cov_S, self$model)


    },


    #' @description
    #' Perform asymptotic statistical inference for ML estimates
    #'
    #' @details
    #' Computes standard errors using the inverse of the information matrix.
    #' Provides z-statistics and p-values based on asymptotic normality.
    #'
    #' @return Invisible self (for method chaining)
    ml_infer = function(){
      theta_ml <- self$estimate$theta
      S <- self$data$cov_S
      N <- self$data$n_row

      ml_infer_estimate <- mlSEM_infer(theta_ml, S, self$model, N, self$estimate)


      self$infer_estimate <- ml_infer_estimate$estimate
      self$infer_estimate$VCOV <- ml_infer_estimate$VCOV
      self$infer_estimate$vcov_effect <- ml_infer_estimate$vcov_effect
    },






    #' @description
    #' Calculate goodness-of-fit statistics
    #'
    #' @param B Integer number of bootstrap replications for Bollen-Stine test (default: 1000).
    #'   Only used when estimator is "svd".
    #'
    #' @details
    #' Computes multiple fit indices including:
    #' \itemize{
    #'   \item Reliability coefficients for reflective blocks (Dillon)
    #'   \item Chi-square test statistic
    #'   \item CFI (Comparative Fit Index)
    #'   \item TLI (Tucker-Lewis Index)
    #'   \item RMSEA (Root Mean Square Error of Approximation)
    #'   \item SRMR (Standardized Root Mean Square Residual)
    #'   \item Information criteria (AIC, BIC, SABIC)
    #'   \item Bollen-Stine bootstrap p-value (for SVD only)
    #' }
    #'
    #' @return Invisible self (for method chaining)

    get_gof = function(B = 1000){

      estimator <- self$estimator

      res_gof <- list()

      # reliability only for relflective block (Dillon)
      if (sum(self$model$mode == "reflective") > 0){
        res_reliability <- reliability('Dillon', self$estimate$lambda, self$estimate$residual_variance)
        res_gof$reliability <- res_reliability
      }


      if (estimator == 'svd' && is.null(self$gof$bollen_stine)){
        bollen_stine <- svdSEM_gof(self$estimate, B)
        res_gof$bollen_stine <- bollen_stine
      } else if (estimator == 'ml'){
        p <- self$model$p
        q <- self$model$q
        r <- self$model$r
        F <- self$gof$F
        N <- self$data$n_row
        S <- self$data$cov_S
        Sigma <- self$estimate$SIGMA_IMPLIED

        chi2 <- chi2sem(p, q, r, F, N)
        res_gof$chi2 <- chi2

        # basline test
        S_baseline <- diag(diag(S))
        F_baseline <- log(det(S_baseline)) + sum(diag(S%*%solve(S_baseline))) - log(det(S)) - NCOL(S)
        baseline <- chi2sem(p, p, 0, F_baseline, N)
        res_gof$baseline <- baseline


        cfi <- 1 - (max(chi2$test - chi2$df, 0)) /
             (max(baseline$test - baseline$df, chi2$test - chi2$df, 0))

        tli <- ( (baseline$test /  baseline$df) - (chi2$test / chi2$df) ) /
               ( (baseline$test /  baseline$df) - 1 )

        res_gof$cfi <- cfi
        res_gof$tli <- tli


        # rmsea
        RMSEA <- rmseasem(chi2$test, chi2$df, N)
        res_gof$RMSEA <- RMSEA
        SRMR <- srmrsem(S, Sigma)
        res_gof$SRMR <- SRMR

        # loglik
        loglik_H0 <- -(N/2)*(p*log(2*pi) + log(det(Sigma)) + sum(diag(solve(Sigma) %*% S)))
        res_gof$info_criteria$loglik_H0 <- loglik_H0
        loglik_H1 <- -(N/2)*(p*log(2*pi) + log(det(S)) + sum(diag(solve(S) %*% S)))
        res_gof$info_criteria$loglik_H1 <- loglik_H1
        AIC <- -2 * loglik_H0 + 2 * q
        res_gof$info_criteria$AIC <- AIC
        BIC <- -2 * loglik_H0 + q * log(N)
        res_gof$info_criteria$BIC <- BIC
        SABIC <- -2 * loglik_H0 + q * log((N + 2) / 24)
        res_gof$info_criteria$SABIC <- SABIC

      }
      self$gof <- c(self$gof, res_gof)

    },



     #' @description
    #' Fit the complete model with inference and goodness-of-fit
    #'
    #' @param infer Logical indicating whether to perform statistical inference (default: FALSE)
    #' @param B Integer number of bootstrap replications for svd (default: 1000)
    #' @param initialisation_svd Logical indicating whether to use SVD initialization
    #'   for ML estimation (default: TRUE). Ignored when estimator is "svd".
    #' @param tol Numeric tolerance for convergence in ML optimization (default: 1e-8)
    #'
    #' @details
    #' This is the main wrapper function that performs:
    #' \enumerate{
    #'   \item Parameter estimation
    #'   \item Statistical inference
    #'   \item Goodness-of-fit assessment
    #' }
    #'
    #' @return Invisible self (for method chaining)
    #'
    #' @examples
    #' \dontrun{
    #' model$fit(estimator = "svd", B = 1000)
    #' model$fit(estimator = "ml", B = 500, initialisation_svd = TRUE)
    #' }

    fit = function(infer = FALSE, B = 1000, initialisation_svd = TRUE, tol = 1e-8){
      estimator <- self$estimator
      self$boot_rep <- B
      if (estimator == 'svd'){
        self$fit_svd()

      } else if(estimator == 'ml'){
        self$fit_ml(initialisation_svd, tol)
        self$get_gof()
      }
      if (infer){
        if (estimator == 'svd'){
          boot_out <- bootstrap_svd(self$estimate, B, verbose = TRUE)
          self$infer_estimate <- boot_out$infer
          self$gof$bollen_stine <- boot_out$gof
          self$get_gof()
        }
      else if(estimator == 'ml'){
          self$ml_infer()
        }
      }



    },

    #' @description
    #' Print comprehensive summary of model estimation results
    #'
    #' @param standardized Logical indicating whether to display standardized estimates (default: FALSE)
    #' @param effect Logical indicating whether to display total and indirect effects (default: FALSE)
    #' @param all_measures Logical indicating whether to display all goodness-of-fit measures (default: FALSE)
    #'
    #' @details
    #' Displays:
    #' \itemize{
    #'   \item Model information (estimator, sample size, number of parameters)
    #'   \item Chi-square test results
    #'   \item Baseline model comparison
    #'   \item Fit indices (CFI, TLI)
    #'   \item Information criteria (AIC, BIC, SABIC)
    #'   \item RMSEA with confidence intervals
    #'   \item SRMR
    #'   \item Bollen-Stine bootstrap results (SVD only)
    #'   \item Parameter estimates with standard errors and p-values
    #' }
    #'
    #' @return Invisible NULL

    summary = function(standardized = FALSE, effect = FALSE, all_measures  = F){

      estimator <- self$estimator



      cat("\n")
      cat(sprintf("%-45s%15s\n", "Estimator", toupper(estimator)))
      cat(sprintf("%-45s%15d\n", "Number of model parameters", sum(self$model$lengths_theta)))
      cat(sprintf("%-45s%15d\n", "Number of observations", self$data$n_row))
      cat(sprintf("%-45s%15d\n", "Degrees of freedom", self$model$dof))
      cat(sprintf("%-45s%15.3f\n", "F", self$gof$F))
      cat(sprintf("%-45s%15.3f\n", "d_LS", self$estimate$T_LS))
      cat("\n")

      if (all_measures){

        if (estimator == 'ml'){

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




          # RMSEA and srmr
          rmsea_val <- self$gof$RMSEA$estimate
          rmsea_ci_lower <- self$gof$RMSEA$CI_lower
          rmsea_ci_upper <- self$gof$RMSEA$CI_upper
          p_rmsea_le_005 <- self$gof$RMSEA$p_close_fit
          p_rmsea_ge_008 <- self$gof$RMSEA$p_notclose_fit
          srmr_val <- self$gof$SRMR
          cat("Model Test User Model :\n\n")
          cat(sprintf("  %-40s%12.3f\n", "Test statistic", testchi2))
          cat(sprintf("  %-40s%12d\n", "Degrees of freedom", dfchi2))
          cat(sprintf("  %-40s%12.3f\n", "P-value (Chi-square)", pvalchi2))
          cat("\n")

          cat("Model Test Baseline Model :\n\n")
          cat(sprintf("  %-40s%12.3f\n", "Test statistic", testbaseline))
          cat(sprintf("  %-40s%12d\n", "Degrees of freedom", dfbaseline))
          cat(sprintf("  %-40s%12.3f\n", "P-value", pvalbaseline))
          cat("\n")

          cat("User Model versus Baseline Model:\n\n")
          cat(sprintf("  %-40s%12.3f\n", "Comparative Fit Index (CFI)", cfi))
          cat(sprintf("  %-40s%12.3f\n", "Tucker-Lewis Index (TLI)", tli))
          cat("\n")

          cat("Root Mean Square Error of Approximation:\n\n")
          cat(sprintf("  %-40s%12.3f\n", "RMSEA", rmsea_val))
          cat(sprintf("  %-40s%12.3f\n", "90 Percent confidence interval - lower", rmsea_ci_lower))
          cat(sprintf("  %-40s%12.3f\n", "90 Percent confidence interval - upper", rmsea_ci_upper))
          cat(sprintf("  %-40s%12.3f\n", "P-value H_0: RMSEA <= 0.050", p_rmsea_le_005))
          cat(sprintf("  %-40s%12.3f\n", "P-value H_0: RMSEA >= 0.080", p_rmsea_ge_008))
          cat("\n")

          cat("Standardized Root Mean Square Residual:\n\n")
          cat(sprintf("  %-40s%12.3f\n", "SRMR", srmr_val))
          cat("\n")

          # Information criteria

          loglik_H0 <-  self$gof$info_criteria$loglik_H0
          loglik_H1 <-  self$gof$info_criteria$loglik_H1
          AIC <-  self$gof$info_criteria$AIC
          BIC <-  self$gof$info_criteria$BIC
          SABIC <-  self$gof$info_criteria$SABIC

          cat("Loglikelihood and Information Criteria:\n\n")
          cat(sprintf("  %-40s%12.3f\n", "Loglikelihood user model (H0)", loglik_H0))
          cat(sprintf("  %-40s%12.3f\n", "Loglikelihood unrestricted model (H1)", loglik_H1))
          cat("\n")
          cat(sprintf("  %-40s%12.3f\n", "Akaike (AIC)", AIC))
          cat(sprintf("  %-40s%12.3f\n", "Bayesian (BIC)", BIC))
          cat(sprintf("  %-40s%12.3f\n", "Sample-size adjusted BIC (SABIC)", SABIC))
          cat("\n")
        }


        if (estimator == 'svd' && !is.null(self$gof$bollen_stine)){
          B <- self$boot_rep
          pvalbs <- self$gof$bollen_stine$pval
          cat("Bootstrap Test (Bollen Stine):\n\n")
          cat(sprintf("  %-40s%12d\n", "Number of bootstrap replications", B))
          cat(sprintf("  %-40s%12.3f\n", "Bollen Stine bootstrap p-value", pvalbs))
          cat("\n")
        }
      }


      # reliability

      if (!is.null(self$gof$reliability)){
        cat("Reliability Coefficients (Dillon):\n\n")
        print(round(self$gof$reliability, 3))
        cat("\n")
      }

      # R2

      if (!is.null(self$estimate$R2)){
          cat("R2:\n\n")
          print(round(self$estimate$R2, 3))
          cat("\n")
      }


      # estimation
      estimate <- formatting_estimate(self$estimate)
      lambda <- estimate$lambda
      if (standardized){
        lambda <- estimate$std_lambda
      }
      beta <- estimate$beta
      gamma <- estimate$gamma
      residualvariance <- estimate$residual_variance
      total_effects <- estimate$total_effects
      indirect_effects <- estimate$indirect_effects
      omega <- estimate$omega




      if (!is.null(self$infer_estimate)){

        # inference estimation
        estimate <- self$infer_estimate
        lambda <- estimate$lambda
        if (standardized){
          lambda <- estimate$std_lambda
        }
        beta<- estimate$beta
        gamma<- estimate$gamma
        residualvariance<- estimate$residual_variance
        total_effects <- estimate$total_effects
        indirect_effects <- estimate$indirect_effects
        if (!is.null(estimate$omega)){
          omega <- estimate$omega
        }


      }





      cat("\nParameter Estimates:\n")
      if (standardized){
        print_link(lambda, 'Standardized Loadings', '=~')
      } else {
        print_link(lambda, 'Loadings', '=~')
      }
      print_link(omega, 'Formative block weight', '<~')
      print_link(rbind(beta, gamma), 'Regression', '~')

      print_variance(residualvariance, 'Residual Variances')

      if (effect){
        print_link(total_effects, 'Total Effects', '~')
        print_link(indirect_effects, 'Indirect Effects', '~')

      }

      # cat("---\nSignif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")


    }

  )
)


