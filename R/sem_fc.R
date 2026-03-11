
#' @import R6
#' @name SemFC
#' @title SemFC Class
#'
#' @description
#' The SEMFC package implements Structural Equation Modeling (SEM) with factors 
#' and composites within the framework of the basic design (Tenenhaus et al, 2025). 
#' 
#' @details
#' The SEMFC package supports the svdSEM and the (restricted) maximum likelihood 
#' estimation methods. svdSEM relies on a non-iterative SVD-based algorithm for 
#' parameter estimation and produces consistent and asymptotically normal 
#' estimators, offering a statistically and computationally sound approach. 
#' 
#' In addtion SEMFC implements the restricted maximum-likelihood (RML-SEM) 
#' approach for the basic design with factors and composites. svdSEM estimates 
#' serve as an initial solution for RML-SEM. The RML-SEM estimator is 
#' implemented using the solnp algorithm (Ye, 1987), available in the Rsolnp package (Ghalanos and Theussl, 2015), 
#' which enables efficient nonlinear optimization with constraints. 
#' 
#' The SEMFC package also provides comprehensive tools for statistical inference, 
#' including bootstrap methods for SVD-based estimation and asymptotic inference 
#' for ML, as well as a wide range of goodness-of-fit measures to evaluate model 
#' fit (chi-square, CFI, TLI, RMSEA, SRMR, AIC, BIC).
#'   
#' \strong{References}
#' 
#' \enumerate{
#' \item Tenenhaus, A., Tenenhaus, M., Dijkstra, T.K. Structural equation modeling 
#' with factors and composites within the framework of the basic design. 
#' Advances in Data Analysis Classification (2025). 
#' \url{https://doi.org/10.1007/s11634-025-00647-4}
#' \item Ye Y (1987) Interior algorithms for linear, quadratic, and 
#' linearly constrained non-linear programming. \href{https://web.stanford.edu/~yyye/YinyuYePhD.pdf}{PhD thesis}, Department of 
#' ESS, Stanford University
#' \item Ghalanos A, Theussl S (2015) Rsolnp: general non-linear 
#' optimization using augmented Lagrange multiplier method. 
#' R package version 1.16. \cr
#' \url{https://CRAN.R-project.org/package=Rsolnp}
#' }
#' 
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
#' @field estimate List containing all estimated model parameters:
#'   \itemize{
#'     \item \code{lambda}: List of loading vectors for each block
#'     \item \code{omega}: List of composite weight vectors for formative blocks
#'     \item \code{beta}: Matrix of structural paths between endogenous variables
#'     \item \code{gamma}: Matrix of structural paths from exogenous to endogenous variables
#'     \item \code{residual_variance}: List of residual variances for each observed variable
#'     \item \code{P_EXO}: Correlation matrix of exogenous latent variables
#'     \item \code{P_ENDO}: Correlation matrix of endogenous latent variables
#'     \item \code{P_IMPLIED}: Implied correlation matrix of all latent variables
#'     \item \code{SIGMA_IMPLIED}: Implied covariance matrix of observed variables
#'     \item \code{std_lambda}: List of standardized loadings for each block
#'     \item \code{std_omega}: List of standardized composite weights for formative blocks
#'     \item \code{psi}: Residual covariance matrix of latent variables
#'     \item \code{R2}: Named vector of R-squared values for endogenous latent variables
#'     \item \code{T_LS}: Least squares fit function value
#'     \item \code{theta}: Numeric vector of all free parameters
#'     \item \code{effect}: List of total and indirect effects between latent variables
#'     \item \code{Ptilde}: First-step correlation matrix estimate (SVD only)
#'   }
#' @field infer_estimate List containing inference results for all parameters:
#'   \itemize{
#'     \item \code{lambda}: Data frame of loadings with SE, z-values, p-values and CI
#'     \item \code{omega}: Data frame of composite weights with SE, z-values, p-values and CI
#'     \item \code{beta}: Data frame of structural paths (endo) with SE, z-values, p-values and CI
#'     \item \code{gamma}: Data frame of structural paths (exo) with SE, z-values, p-values and CI
#'     \item \code{residual_variance}: Data frame of residual variances with SE, z-values, p-values and CI
#'     \item \code{total_effects}: Data frame of total effects with SE, z-values, p-values and CI
#'     \item \code{indirect_effects}: Data frame of indirect effects with SE, z-values, p-values and CI
#'     \item \code{VCOV}: Variance-covariance matrix of parameter estimates (ML only)
#'     \item \code{vcov_effect}: Variance-covariance matrix of effect estimates (ML only)
#'   }
#' @field boot_rep Integer number of bootstrap samples
#' @field gof List containing goodness-of-fit statistics:
#'   \itemize{
#'     \item \code{F}: Value of the fit function at the optimal solution
#'     \item \code{reliability}: Named vector of reliability coefficients (Dillon-Goldstein rho) for reflective blocks
#'     \item \code{chi2}: Chi-square test statistic
#'     \item \code{df}: Degrees of freedom for chi-square test
#'     \item \code{pvalue}: P-value of the chi-square test
#'     \item \code{CFI}: Comparative Fit Index
#'     \item \code{TLI}: Tucker-Lewis Index
#'     \item \code{RMSEA}: Root Mean Square Error of Approximation
#'     \item \code{RMSEA_CI}: 90\% confidence interval for RMSEA
#'     \item \code{SRMR}: Standardized Root Mean Square Residual
#'     \item \code{AIC}: Akaike Information Criterion
#'     \item \code{BIC}: Bayesian Information Criterion
#'     \item \code{SABIC}: Sample-size Adjusted BIC
#'     \item \code{bollen_stine}: Bollen-Stine bootstrap p-value (SVD only)
#'   }
#'
#' @examples
#' data("ECSI")
#' ECSI = ECSI/10
#' A = list(IMAG = ECSI[, 1:5],
#'         CUEX = ECSI[, 6:8],
#'         PERQ = ECSI[, 9:15],
#'         PERV = ECSI[, 16:17],
#'         CUSA = ECSI[, 18:20],
#'         CUSCO = ECSI[, 21, drop = FALSE],
#'         CUSL = ECSI[, 22:24])
#'
#' C <- matrix(c(0, 1, 0, 0, 1, 0, 1,
#'               0, 0, 1, 1, 1, 0, 0,
#'               0, 0, 0, 1, 1, 0, 0,
#'               0, 0, 0, 0, 1, 0, 0,
#'               0, 0, 0, 0, 0, 1, 1,
#'               0, 0, 0, 0, 0, 0, 1,
#'               0, 0, 0, 0, 0, 0, 0), 7, 7, byrow = TRUE)
#'
#' colnames(C) = rownames(C) = names(A)
#'
#' mode = rep("reflective", 7) ; mode[6] = "formative"
#'
#' sem_model <- SemFC$new(data = A,
#' relation_matrix = C,
#' mode = mode,
#' scale = FALSE,
#' estimator = "svd")
#'
#' sem_model$fit(infer = TRUE, B = 100)
#'
#' sem_model$summary(standardized = TRUE, effect = TRUE, all_measures = TRUE)
#'
#' estimates <- sem_model$parameterEstimates(standardized = TRUE)
#'
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
#' Initialize SemFC object with data and model specification
#'
#' @param data A list that contains J blocks of indicator variables. 
#'   Blocks are either reflective or formative. Each block should be a data 
#'   frame or matrix with rows as observations and columns as indicators.
#' @param relation_matrix Square connection matrix (J x J) defining 
#'   structural connection between latent variables (1 = connection exists, 
#'   0 = no connection).
#' @param mode Character vector of length J specifying the type of 
#'   measurement model ("reflective" or "formative") for each block. 
#'   (default: rep("reflective", J))
#' @param estimator Character string specifying the estimation method: 
#'   "svd" or "ml" (default: "ml").
#' @param scale Logical indicating whether to standardize the input data or 
#'   not (default: FALSE)
#' @param bias Logical indicating whether to apply bias correction in 
#'   covariance estimation (default: FALSE)
#'
#' @examples
#' data(ECSI)
#' ECSI = ECSI/10
#' A = list(CUSTOMER_E = ECSI[, c("CUEX1", "CUEX2", "CUEX3")],
#'          PERC_QUAL  = ECSI[, c("PERQ1", "PERQ2", "PERQ3", 
#'                                "PERQ4", "PERQ5", "PERQ6", "PERQ7")],
#'          PERC_VALUE = ECSI[, c("PERV1", "PERV2")],
#'          CUSTOMER_S = ECSI[, c("CUSA1", "CUSA2", "CUSA3")],
#'          CUSTOMER_L = ECSI[, c("CUSL1", "CUSL2", "CUSL3")])
#'
#' C = matrix(c(0, 0, 0, 0, 0,
#'              1, 0, 0, 0, 0,
#'              1, 1, 0, 0, 0,
#'              1, 1, 1, 0, 0,
#'              0, 0, 0, 1, 0), 5, 5, 
#'              byrow = FALSE)
#'              
#' colnames(C) = rownames(C) = names(A)
#'
#' sem_model <- SemFC$new(data = A,
#' relation_matrix = C,
#' mode = rep("reflective", 5),
#' scale = FALSE,
#' estimator = "svd")
#'
#' @return A new `SemFC` object

initialize = function(data, relation_matrix, mode, estimator = 'ml', 
                      scale = FALSE, bias = FALSE) {

      self$estimator <- estimator
      init <- get_parameter_model_sem(data, mode, relation_matrix, bias)
      self$model <- init$model
      self$data <- init$data

      self$model$scale <- scale
      self$model$bias <- bias

    },

#' @description
#' Fit the full model with inference and goodness-of-fit.
#'
#' @param infer Logical indicating whether to perform statistical inference 
#'   or not (default: FALSE)
#' @param B Integer number of bootstrap samples for svd (default: 1000)
#' @param initialization Character string or numeric vector specifying the 
#'   initialization method for ML. This argument is ignored when estimator 
#'   is "svd".
#'   \itemize{
#'     \item \code{"svd"} (default): use svdSEM estimate as starting values
#'     \item \code{"random"}: use random starting values
#'     \item \code{numeric vector}: use the provided vector as starting 
#'        values. The length of this vector equals the total number of model 
#'        parameters.
#'   }
#' @param tol Numeric tolerance for convergence in ML optimization 
#'   (default: 1e-8)
#'
#' @details
#' This is the main wrapper function that performs:
#' \enumerate{
#'   \item Parameter estimation for svdSEM ("svd") or ML ("ml").
#'   \item Statistical inference (bootstrap for svdSEM, asymptotic for ML)
#'   \item Goodness-of-fit assessment (chi-square, CFI, TLI, RMSEA, SRMR, 
#'     AIC, BIC)
#' }
#' 
#' @examples
#' data(ECSI)
#' ECSI = ECSI/10
#' A = list(CUSTOMER_E = ECSI[, c("CUEX1", "CUEX2", "CUEX3")],
#'          PERC_QUAL  = ECSI[, c("PERQ1", "PERQ2", "PERQ3", 
#'                                "PERQ4", "PERQ5", "PERQ6", "PERQ7")],
#'          PERC_VALUE = ECSI[, c("PERV1", "PERV2")],
#'          CUSTOMER_S = ECSI[, c("CUSA1", "CUSA2", "CUSA3")],
#'          CUSTOMER_L = ECSI[, c("CUSL1", "CUSL2", "CUSL3")])
#'
#' C = matrix(c(0, 0, 0, 0, 0,
#'              1, 0, 0, 0, 0,
#'              1, 1, 0, 0, 0,
#'              1, 1, 1, 0, 0,
#'              0, 0, 0, 1, 0), 5, 5, 
#'              byrow = FALSE)
#'              
#' colnames(C) = rownames(C) = names(A)
#'
#' sem_model_svd <- SemFC$new(data = A,
#'                            relation_matrix = C,
#'                            mode = rep("reflective", 5), 
#'                            scale = FALSE, 
#'                            estimator = "svd")
#'
#' sem_model_svd$fit(infer = TRUE, B = 100)
#'
#' sem_model_ml <- SemFC$new(data = A, 
#'                           relation_matrix = C, 
#'                           mode = rep("reflective", 5), 
#'                           scale = FALSE, 
#'                           estimator = "ml")
#'
#' sem_model_ml$fit(infer = TRUE, 
#'                  initialization = "svd", 
#'                  tol = 1e-04)
#'

fit = function(infer = FALSE, B = 1000, initialization = 'svd', tol = 1e-8){
      estimator <- self$estimator
      self$boot_rep <- B
      if (estimator == 'svd'){
        private$fit_svd()

      } else if(estimator == 'ml'){
        private$fit_ml(initialization, tol)
        private$get_gof()
      }
      if (infer){
        if (estimator == 'svd'){
          private$svd_infer(B, verbose = TRUE)
          private$get_gof()
        }
      else if(estimator == 'ml'){
          private$ml_infer()
        }
      }
    },

#' @description
#' Print comprehensive summary of model estimation results.
#'
#' @param standardized Logical indicating whether to report standardized 
#'   estimates (default: FALSE)
#' @param effect Logical indicating whether to report total and indirect 
#'   effects (default: FALSE)
#' @param all_measures Logical indicating whether to report all 
#'   goodness-of-fit measures (default: FALSE)
#'
#' @details
#' Elements that are reported in the summary include:
#' \itemize{
#'   \item Model information (estimator, sample size, number of parameters)
#'   \item Chi-square test results
#'   \item Baseline model comparison
#'   \item Fit indices (CFI, TLI)
#'   \item Information criteria (AIC, BIC, SABIC)
#'   \item RMSEA with confidence intervals
#'   \item SRMR
#'   \item Bollen-Stine bootstrap results (for SVD only)
#'   \item Parameter estimates with standard errors and p-values
#' }
#'
#' @examples
#' data(ECSI)
#' data(ECSI)
#' ECSI = ECSI/10
#' A = list(CUSTOMER_E = ECSI[, c("CUEX1", "CUEX2", "CUEX3")],
#'          PERC_QUAL  = ECSI[, c("PERQ1", "PERQ2", "PERQ3", 
#'                                "PERQ4", "PERQ5", "PERQ6", "PERQ7")],
#'          PERC_VALUE = ECSI[, c("PERV1", "PERV2")],
#'          CUSTOMER_S = ECSI[, c("CUSA1", "CUSA2", "CUSA3")],
#'          CUSTOMER_L = ECSI[, c("CUSL1", "CUSL2", "CUSL3")])
#'
#' C = matrix(c(0, 0, 0, 0, 0,
#'              1, 0, 0, 0, 0,
#'              1, 1, 0, 0, 0,
#'              1, 1, 1, 0, 0,
#'              0, 0, 0, 1, 0), 5, 5, 
#'              byrow = FALSE)
#'              
#' colnames(C) = rownames(C) = names(A)
#'
#' sem_model <- SemFC$new(data = A, 
#'                        relation_matrix = C, 
#'                        mode = rep("reflective", 5), 
#'                        scale = FALSE, 
#'                        estimator = "svd")
#'
#' sem_model$fit(infer = TRUE, B = 100)
#'
#' sem_model$summary(standardized = TRUE, 
#'                   effect = TRUE, 
#'                   all_measures = TRUE)
#'                   

    summary = function(standardized = F, effect = FALSE, all_measures  = F){

      estimator <- self$estimator
      print_model(estimator, sum(self$model$lengths_theta), self$data$n_row,
                  self$model$dof, self$gof$F, self$estimate$T_LS)

      print_gof(all_measures, estimator, self$gof, self$boot_rep, self$estimate$R2)

      # estimation
      estimate <- formatting_estimate(self$estimate)
      if (!is.null(self$infer_estimate)){
        estimate <- self$infer_estimate
      }
      print_estimates(estimate, standardized, effect)

    },

#' @description 
#' Extract parameter estimates from the fitted SemFC model
#'
#' @param standardized Logical indicating whether to include standardized
#'   estimates (default: FALSE). When TRUE, adds a `std.all` column with
#'   fully standardized coefficients.
#'
#' @details
#' Returns a data frame containing all estimated parameters including:
#' \itemize{
#'   \item \code{lambda}: Loading vectors for reflective and formative 
#'     blocks
#'   \item \code{omega}: Composite weights for formative blocks only
#'   \item \code{beta}: Structural coefficients between endogenous 
#'     latent variables
#'   \item \code{gamma}: Structural coefficients from exogenous to 
#'     endogenous latent variables
#'   \item \code{residualvariance}: Residual variances for observed 
#'     variables in reflective blocks
#' }
#'
#' If statistical inference has been performed, the returned estimates
#' include standard errors, z-values, p-values, and confidence intervals.
#'
#' @return A data frame with all parameter estimates. Columns include:
#' \itemize{
#'   \item \code{lhs}: Left-hand side variable
#'   \item \code{op}: Operator
#'   \item \code{rhs}: Right-hand side variable
#'   \item \code{est}: Point estimate
#'   \item \code{se}: Standard error (if inference was performed)
#'   \item \code{z}: Z-statistic (if inference was performed)
#'   \item \code{pvalue}: P-value (if inference was performed)
#'   \item \code{std.all}: Standardized estimate (if standardized = TRUE)
#' }
#'
#' @examples
#' data(ECSI)
#' ECSI = ECSI/10
#' A = list(CUSTOMER_E = ECSI[, c("CUEX1", "CUEX2", "CUEX3")],
#'          PERC_QUAL  = ECSI[, c("PERQ1", "PERQ2", "PERQ3", 
#'                                "PERQ4", "PERQ5", "PERQ6", "PERQ7")],
#'          PERC_VALUE = ECSI[, c("PERV1", "PERV2")],
#'          CUSTOMER_S = ECSI[, c("CUSA1", "CUSA2", "CUSA3")],
#'          CUSTOMER_L = ECSI[, c("CUSL1", "CUSL2", "CUSL3")])
#'
#' C = matrix(c(0, 0, 0, 0, 0,
#'              1, 0, 0, 0, 0,
#'              1, 1, 0, 0, 0,
#'              1, 1, 1, 0, 0,
#'              0, 0, 0, 1, 0), 5, 5, 
#'              byrow = FALSE)
#'              
#' colnames(C) = rownames(C) = names(A)
#'
#' sem_model <- SemFC$new(data = A, 
#'                        relation_matrix = C,
#'                        mode = rep("reflective", 5), 
#'                        scale = FALSE, 
#'                        estimator = "svd")
#'
#' sem_model$fit(infer = TRUE, B = 100)
#'
#' estimates <- sem_model$parameterEstimates(standardized = TRUE)
#'

    parameterEstimates = function(standardized = FALSE){
      estimate <- formatting_estimate(self$estimate)
      if (!is.null(self$infer_estimate)){
        estimate <- self$infer_estimate
      }
      lambda <- estimate$lambda
      residualvariance <- estimate$residual_variance
      beta <- estimate$beta
      gamma <- estimate$gamma
      omega <- estimate$omega


      table_estimate <- rbind(lambda, omega, beta, gamma, residualvariance)
      rownames(table_estimate) <- NULL
      if (!standardized){
        table_estimate$std.all <- NULL
      }

      return(table_estimate)

    },

#' @description
#' Check for improper solutions in the fitted SemFC model
#'
#' @details
#' Detects improper solutions for:
#' \itemize{
#'   \item \code{reliability_coef}: reliability coefficients for each block
#'   \item \code{p_implied}: implied correlation matrix between latent 
#'     variables
#'   \item \code{residual_variance}: residual variances
#'   \item \code{std_lambda}: standardized loadings
#'   \item \code{sigma_implied}: implied covariance matrix of observed 
#'     variables
#'   \item \code{r2}: r-squared values for endogenous latent variables
#'   \item \code{psi}: residual covariance matrix of latent variables
#'   \item \code{p_tilde}: intermediate estimation of correlation matrix 
#'     between latent variables (only for SVDSEM)
#' }
#'
#' @return logical vector indicating improper solutions for each of the 
#'   following criteria:
#' \describe{
#'   \item{reliability_coef}{\code{TRUE} if any reliability coefficient is 
#'     outside (0, 1).}
#'   \item{rho_jh}{\code{TRUE} if any correlation in P_IMPLIED is outside 
#'     (-1, 1).}
#'   \item{p_implied}{\code{TRUE} if p_implied is not positive definite.}
#'   \item{theta_jh}{\code{TRUE} if any residual variance is negative.}
#'   \item{std_lambda}{\code{TRUE} if any standardized loading is outside 
#'     (-1, 1).}
#'   \item{sigma_implied}{\code{TRUE} if sigma_implied is not positive 
#'     definite.}
#'   \item{r2}{\code{TRUE} if any r-squared is outside (0, 1).}
#'   \item{psi}{\code{TRUE} if psi is not positive definite.}
#'   \item{p_tilde}{\code{TRUE} if p_tilde is not positive definite (only 
#'     for svdSEM).}
#' }
#'
#' @examples
#' data(ECSI)
#' ECSI = ECSI/10
#' A = list(CUSTOMER_E = ECSI[, c("CUEX1", "CUEX2", "CUEX3")],
#'          PERC_QUAL  = ECSI[, c("PERQ1", "PERQ2", "PERQ3", 
#'                                "PERQ4", "PERQ5", "PERQ6", "PERQ7")],
#'          PERC_VALUE = ECSI[, c("PERV1", "PERV2")],
#'          CUSTOMER_S = ECSI[, c("CUSA1", "CUSA2", "CUSA3")],
#'          CUSTOMER_L = ECSI[, c("CUSL1", "CUSL2", "CUSL3")])
#'
#' C = matrix(c(0, 0, 0, 0, 0,
#'              1, 0, 0, 0, 0,
#'              1, 1, 0, 0, 0,
#'              1, 1, 1, 0, 0,
#'              0, 0, 0, 1, 0), 5, 5, 
#'              byrow = FALSE)
#'              
#' colnames(C) = rownames(C) = names(A)
#'
#' sem_model <- SemFC$new(data = A, 
#'                        relation_matrix = C, 
#'                        mode = rep("reflective", 5), 
#'                        scale = FALSE, 
#'                        estimator = "svd")
#'
#' sem_model$fit(infer = TRUE, B = 100)
#'
#' improper_results <- sem_model$check_improper()
#'

    check_improper = function(){
      if (is.null(self$estimate) || length(self$estimate) == 0L) {
        stop("Model has not been fitted yet. Call `fit()` before `check_improper()`.", call. = FALSE)
      }
      return(improper(self$estimate))
    }
  ),

  private = list(
#
# ' @description
# ' Fit the model using Singular Value Decomposition (SVD) method
# '
# ' @details
# ' Estimates model parameters using SVD-based approach which is computationally
# ' efficient and provides good starting values for ML estimation.
# '
# ' @return Invisible self (for method chaining)
# ' @keywords internal

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
    
# ' @description
# ' Perform statistical inference for SVD estimates using bootstrap
# '
# ' @param B Integer number of bootstrap samples (default: 1000)
# ' @param verbose Logical indicating whether to print progress messages (default: TRUE)
# '
# ' @details
# ' Uses non-parametric bootstrap to estimate standard errors, confidence intervals,
# ' and p-values for all model parameters.
# '
# ' @return Invisible self (for method chaining)
# ' @keywords internal

    svd_infer = function(B = 1000, verbose = TRUE){

      boot_out <- svdsem_infer(self$estimate, B, verbose = verbose)
      self$infer_estimate <- boot_out$result$infer
      self$gof$bollen_stine <- boot_out$gof


    },

# ' @description
# ' Fit the model using Maximum Likelihood (ML) estimation
# '
# ' @param initialization Character string or numeric vector specifying 
# '    initialization method:
# '   \itemize{
# '     \item \code{"svd"} (default): Use svdSEM estimate as starting values
# '     \item \code{"random"}: Use random starting values
# '     \item \code{numeric vector}: use the provided vector as starting 
# '              values. This vector must have length equal to the total 
# '              number of model parameters
# '   }
# ' @param tol Numeric tolerance for convergence in optimization (default: 1e-8)
# '
# ' @details
# ' Uses numerical optimization (via SOLNP) to minimize the ML fit function.
# ' SVD initialization is recommended for better convergence.
# '
# ' @return Invisible self (for method chaining)
# ' @keywords internal

  fit_ml = function(initialization = 'svd', tol) {

      len_theta <- sum(self$model$lengths_theta)

      initial_params <- if (is.numeric(initialization)) {
        if (length(initialization) != len_theta) {
          stop("Length of provided initial parameters does not match the number of model parameters. ",
               "Expected length: ", len_theta)
        }
        initialization
      } else {
        initialization <- match.arg(initialization, c('svd', 'random'))
        switch(initialization,
          svd = {
            private$fit_svd()
            self$estimator <- 'ml'
            self$estimate$theta
          },
          random = runif(len_theta)
        )
      }

      ml_sol <- mlSEM(initial_params, self$data$cov_S, self$model, tol)
      theta_ml <- ml_sol$pars
      self$estimate <- lvm_ml(x = theta_ml, model = self$model, jac = F)
      self$estimate$T_LS <- d_LS(self$data$cov_S, self$estimate$SIGMA_IMPLIED)

      var_MVs <- lapply(self$data$data, 
                        function(x) 
                          diag(cov2(x, bias = self$model$bias)))
      std_lambda <- mapply("/", self$estimate$lambda, lapply(var_MVs, sqrt),  
                           SIMPLIFY = FALSE)
      self$estimate$std_lambda <- std_lambda
      std_omega <- mapply(function(Sjj, lambda_j) solve(Sjj) %*% lambda_j,
                         self$data$S_diag_composites, 
                         std_lambda[self$model$mode == "formative"],
                         SIMPLIFY = FALSE)
      names(std_omega) <- names(std_lambda[self$model$mode == "formative"])
      self$estimate$std_omega <- std_omega
      self$estimate$theta <- theta_ml
      self$estimate$effect <- compute_effect(self$estimate$beta, 
                                             self$estimate$gamma)
      self$gof$F <- F1(theta_ml, self$data$cov_S, self$model)
    },


# ' @description
# ' Perform asymptotic statistical inference for restricted Maximum Likelihood estimates
# '
# ' @details
# ' Computes standard errors using the inverse of the information matrix.
# ' Provides z-statistics and p-values based on asymptotic normality.
# '
# ' @return Invisible self (for method chaining)
# ' @keywords internal

    ml_infer = function(){
      theta_ml <- self$estimate$theta
      S <- self$data$cov_S
      N <- self$data$n_row

      ml_infer_estimate <- mlSEM_infer(theta_ml, S, self$model, N, self$estimate)


      self$infer_estimate <- ml_infer_estimate$estimate
      self$infer_estimate$VCOV <- ml_infer_estimate$VCOV
      self$infer_estimate$vcov_effect <- ml_infer_estimate$vcov_effect
    },


# ' @description
# ' Calculate goodness-of-fit statistics
# '
# ' @param B Integer number of bootstrap samples for Bollen-Stine test (default: 1000).
# '   Only used when estimator is "svd".
# '
# ' @details
# ' Computes multiple fit indices including:
# ' \itemize{
# '   \item Reliability coefficients for reflective blocks (Dillon)
# '   \item Chi-square test statistic
# '   \item CFI (Comparative Fit Index)
# '   \item TLI (Tucker-Lewis Index)
# '   \item RMSEA (Root Mean Square Error of Approximation)
# '   \item SRMR (Standardized Root Mean Square Residual)
# '   \item Information criteria (AIC, BIC, SABIC)
# '   \item Bollen-Stine bootstrap p-value (for SVD only)
# ' }
# '
# ' @return Invisible self (for method chaining)
# ' @keywords internal

    get_gof = function(B = 1000){

      estimator <- self$estimator

      res_gof <- list()

      # reliability only for relflective block (Dillon)
      if (sum(self$model$mode == "reflective") > 0){
        res_reliability <- reliability('Dillon', self$estimate$lambda, 
                                       self$estimate$residual_variance)
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

        res_gof <- c(res_gof, semML_gof(p, q, r, F, N, S, Sigma))

      }
      self$gof <- c(self$gof, res_gof)

    }
  )
)


