
#' @import R6
#' @name SemFC
#' @title SemFC: Structural Equation Modeling (SEM) with factors and composites
#' within the framework of the basic design
#'
#' @description
#' The SEMFC package implements Structural Equation Modeling (SEM) with factors
#' and composites within the framework of the basic design (Tenenhaus et al, 2025).
#'
#' @details
#' The SEMFC package supports the svdSEM and the (restricted) maximum likelihood
#' estimation methods.
#'
#' \itemize{
#'
#' \item svdSEM relies on a non-iterative SVD-based algorithm for
#' parameter estimation and produces consistent and asymptotically normal
#' estimators, offering a statistically and computationally sound approach.
#'
#' \item The restricted maximum-likelihood (RML-SEM) approach for the basic
#' design with factors and composites is also implemented within SEMFC.
#' svdSEM estimates serve as an initial solution for RML-SEM. The RML-SEM
#' estimator is implemented using the solnp algorithm (Ye, 1987), available
#' in the Rsolnp package (Ghalanos and Theussl, 2015), which enables efficient
#' nonlinear optimization with constraints.}
#'
#' The SEMFC package also provides comprehensive tools for statistical inference,
#' including bootstrap methods for SVD-based estimation and asymptotic inference
#' for ML, as well as a wide range of goodness-of-fit measures to evaluate model
#' fit (chi-square, CFI, TLI, RMSEA, SRMR, AIC, BIC).
#'
#' \strong{References}
#'
#' \enumerate{
#' \item Tenenhaus, A., Tenenhaus, M., Dijkstra, T.K. Structural equation 
#' modeling with factors and composites within the framework of the basic design.
#' Advances in Data Analysis Classification (2025).
#' \url{https://doi.org/10.1007/s11634-025-00647-4}
#' \item Ye Y (1987) Interior algorithms for linear, quadratic, and
#' linearly constrained non-linear programming. 
#' \href{https://web.stanford.edu/~yyye/YinyuYePhD.pdf}{PhD thesis}, Department 
#' of ESS, Stanford University
#' \item Ghalanos A, Theussl S (2015) Rsolnp: general non-linear
#' optimization using augmented Lagrange multiplier method.
#' R package version 1.16. \cr
#' \url{https://CRAN.R-project.org/package=Rsolnp}
#' }
#'
#' @export
SemFC <- R6Class(
  "SemFC",
  public = list(


#' @description
#' Initialize SemFC object with data and model specification
#'
#' @param data A list that contains \code{J} blocks of indicator variables.
#'   Blocks are either reflective or formative. Each block should be a data
#'   frame or matrix with rows as observations and columns as indicators.
#' @param relation_matrix Square 0/1 connection matrix (\code{J} x \code{J})
#'   defining structural connection between latent variables.
#'   \code{relation_matrix[i, j]} = 1 indicates a structural
#'   connection from latent variable \code{i} to latent variable \code{j}; 
#'   and equal 0 otherwise.
#' @param mode Character vector of length \code{J} specifying the type of
#'   measurement model (\code{"reflective"} or \code{"formative"}) for each block.
#'   Blocks specified as
#'   \code{"reflective"} are modeled with latent factors, while blocks specified
#'   as \code{"formative"} are modeled with composites. 
#'   (default: \code{rep("reflective", J)})
#' @param estimator Character string specifying the estimation method:
#'   \code{"svd"} or \code{"ml"} (default: \code{"ml"}).
#' @param scale Logical indicating whether to standardize the input data or
#'   not. When \code{TRUE}, all variables are standardized to zero-mean and
#'   unit variance before estimation. When \code{FALSE}, centered data is used.
#'   (default: \code{FALSE})
#' @param bias Logical indicating whether to apply bias correction in
#'   covariance estimation. When \code{TRUE}, covariance matrix is computed
#'   with division by \code{n} (biased estimator);  when \code{FALSE}, division
#'   by \code{n-1} is used (unbiased estimator). (default: \code{FALSE})
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
#' # Initialize SemFC object with data, model specification 
#' # and svd-based estimation method.
#' 
#' sem_model <- SemFC$new(data = A, 
#'                        relation_matrix = C, 
#'                        mode = rep("reflective", 5), 
#'                        scale = FALSE, estimator = "svd")
#'
#' @return A new `SemFC` object

initialize = function(data, relation_matrix, 
                      mode = rep('reflective', length(data)), 
                      estimator = 'ml',
                      scale = FALSE, bias = FALSE) {

      private$.estimator <- estimator
      init <- get_parameter_model_sem(data, mode, relation_matrix, bias)
      private$.model <- init$model
      private$.data <- init$data

      private$.model$scale <- scale
      private$.model$bias <- bias

    },

#' @description
#' Fit the full model with inference and goodness-of-fit.
#'
#' @param infer Logical indicating whether to perform statistical inference
#'   or not. When TRUE, bootstrap inference is performed for svdSEM and
#'   asymptotic inference is performed for ML. (default: \code{FALSE})
#' @param B Integer number of bootstrap samples for svd (default: 1000)
#' @param initialization Character string or numeric vector specifying the
#'   initialization strategy for the ML algorithm. This argument is ignored 
#'   when estimator is \code{"svd"}. 
#'   \itemize{
#'     \item \code{"svd"} (default): use svdSEM estimate as starting values
#'     \item \code{"random"}: use random starting values
#'     \item \code{numeric vector}: use the provided vector as starting
#'        values. The length of this vector equals the total number of model
#'        parameters.
#'   }
#' @param tol Numeric tolerance for convergence in ML optimization
#'   (default: \code{1e-8})
#'
#' @details
#' This is the main wrapper function that performs:
#' \enumerate{
#'   \item Parameter estimation for svdSEM (\code{"svd"}) or ML (\code{"ml"}).
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
#' # Initialize SemFC object with data, model specification 
#' # and svd-based estimation method.
#' 
#' sem_model_svd <- SemFC$new(data = A,
#'                            relation_matrix = C,
#'                            mode = rep("reflective", 5),
#'                            scale = FALSE,
#'                            estimator = "svd")
#'                            
#'# Fit the full model with bootstrap inference and goodness-of-fit.
#'
#' sem_model_svd$fit(infer = TRUE, B = 100)
#' 
#' # Initialize SemFC object with data model specification 
#' # and ml estimation method.
#' 
#' sem_model_ml <- SemFC$new(data = A,
#'                           relation_matrix = C,
#'                           mode = rep("reflective", 5),
#'                           scale = FALSE,
#'                           estimator = "ml")
#'
#'# Fit the full model with asymptotic inference and goodness-of-fit.
#' sem_model_ml$fit(infer = TRUE,
#'                  initialization = "svd",
#'                  tol = 1e-04)
#'

    fit = function(infer = FALSE, B = 1000, initialization = 'svd', tol = 1e-8){
      estimator <- private$.estimator
      private$.boot_rep <- B
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
#' Print comprehensive summary of the fitted SemFC object.
#'
#' @param standardized Logical indicating whether to report standardized
#'   estimates (default: \code{FALSE}).
#' @param effect Logical indicating whether to report total and indirect
#'   effects (default: \code{FALSE})
#' @param all_measures Logical indicating whether to report all
#'   goodness-of-fit measures. When \code{FALSE}, only a subset of key fit
#'   indices is reported. (default: \code{FALSE})
#'
#' @details
#' Elements that are reported in the summary include:
#' \itemize{
#'   \item Model information (estimator, sample size, number of parameters)
#'   \item Chi-square test: This is the primary test of our model. If the test 
#'   is not significant (p > 0.05), it suggests that the model-implied 
#'   covariance matrix does not significantly differ from the sample 
#'   covariance matrix, indicating good fit. 
#'   \item Baseline model comparison
#'   \item Fit indices: Comparative Fit Index (CFI) and Tucker-Lewis Index 
#'   (TLI). Values above 0.95 indicate good fit.
#'   \item Root Mean Square Error of Approximation (RMSEA) with confidence 
#'   intervals: values below 0.05 indicate good fit.
#'   \item Standardized Root Mean Square Residual (SRMR): values below 0.06 
#'   indicate good fit.
#'   \item Information criteria (AIC, BIC) and Sample-size adjusted BIC (SABIC)
#'   \item Bollen-Stine bootstrap results (for SVD only)
#'   \item Parameter estimates with standard errors and p-values
#' }
#' There is a lot of controversy about the use (and misuse) of these fit 
#' indices. Current practice is to report: chi-square value, degree of freedom 
#' and pvalue, RMSEA, CFI and SRMR.
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
#' # Initialize SemFC object with data, model specification 
#' # and svd-based estimation method.
#' 
#' sem_model <- SemFC$new(data = A,
#'                        relation_matrix = C,
#'                        mode = rep("reflective", 5),
#'                        scale = FALSE,
#'                        estimator = "svd")
#'                        
#'# Fit the full model with bootstrap inference and goodness-of-fit.
#'
#' sem_model$fit(infer = TRUE, B = 100)
#' 
#' # Print comprehensive summary of the fitted SemFC object.
#'
#' sem_model$summary(standardized = TRUE,
#'                   effect = TRUE,
#'                   all_measures = TRUE)
#'

    summary = function(standardized = F, effect = FALSE, all_measures  = F){

      estimator <- private$.estimator
      print_model(estimator, sum(private$.model$lengths_theta), private$.data$n_row,
                  private$.model$dof, private$.gof$F, private$.estimate$T_LS)

      print_gof(all_measures, estimator, private$.gof, private$.boot_rep, private$.estimate$r2)

      # estimation
      estimate <- formatting_estimate(private$.estimate)
      if (!is.null(private$.infer_estimate)){
        estimate <- private$.infer_estimate
      }
      print_estimates(estimate, standardized, effect)

    },

#' @description
#' Extract parameter estimates from the fitted SemFC model
#'
#' @param standardized Logical indicating whether to include standardized
#'   estimates ((default: \code{FALSE}). When \code{TRUE}, adds a `std.all`
#'   column with fully standardized coefficients.
#'
#' @details
#' Returns a data frame containing all estimated parameters including:
#' \itemize{
#'   \item \code{lambda}: Loading vectors for reflective and formative
#'     blocks.
#'   \item \code{omega}: Composite weights for formative blocks only
#'   \item \code{beta}: Structural coefficients between endogenous
#'     latent variables.
#'   \item \code{gamma}: Structural coefficients from exogenous to
#'     endogenous latent variables.
#'   \item \code{residualvariance}: Residual variances for observed
#'     variables. Only applies for reflective blocks, as formative blocks do
#'     not have residual variances.
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
#' # Initialize SemFC object with data, model specification 
#' # and svd-based estimation method.
#'
#' sem_model <- SemFC$new(data = A,
#'                        relation_matrix = C,
#'                        mode = rep("reflective", 5),
#'                        scale = FALSE,
#'                        estimator = "svd")
#'                        
#'# Fit the full model with bootstrap inference and goodness-of-fit.
#'
#' sem_model$fit(infer = TRUE, B = 100)
#' 
#' # Extract parameter estimates from the fitted SemFC model, including 
#' # standardized solutions. 
#'
#' estimates <- sem_model$parameterEstimates(standardized = TRUE)
#'

    parameterEstimates = function(standardized = FALSE){
      estimate <- formatting_estimate(private$.estimate)
      if (!is.null(private$.infer_estimate)){
        estimate <- private$.infer_estimate
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
#'     between latent variables (only for svdSEM)
#' }
#'
#' @return logical vector indicating improper solutions for each of the
#'   following criteria:
#' \describe{
#'   \item{\code{reliability_coef}}{\code{TRUE} if any reliability coefficient is
#'     outside (0, 1).}
#'   \item{\code{rho_jh}}{\code{TRUE} if any correlation in \code{p_implied} is outside
#'     (-1, 1).}
#'   \item{\code{p_implied}}{\code{TRUE} if \code{p_implied} is not positive definite.}
#'   \item{\code{theta_jh}}{\code{TRUE} if any residual variance is negative.}
#'   \item{\code{std_lambda}}{\code{TRUE} if any standardized loading is outside
#'     (-1, 1).}
#'   \item{\code{sigma_implied}}{\code{TRUE} if \code{sigma_implied} is not positive
#'     definite.}
#'   \item{\code{r2}}{\code{TRUE} if any r-squared is outside (0, 1).}
#'   \item{\code{psi}}{\code{TRUE} if \code{psi} is not positive definite.}
#'   \item{\code{p_tilde}}{\code{TRUE} if \code{p_tilde} is not positive definite (only
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
#' # Initialize SemFC object with data, model specification 
#' # and svd-based estimation method.
#'
#' sem_model <- SemFC$new(data = A,
#'                        relation_matrix = C,
#'                        mode = rep("reflective", 5),
#'                        scale = FALSE,
#'                        estimator = "svd")
#'                        
#'# Fit the full model with bootstrap inference and goodness-of-fit.
#'
#' sem_model$fit(infer = TRUE, B = 100)
#' 
#' # Check for improper solutions in the fitted SemFC model.
#'
#' improper_results <- sem_model$check_improper()
#'

    check_improper = function(){
      if (is.null(private$.estimate) || length(private$.estimate) == 0L) {
        stop("Model has not been fitted yet. Call `fit()` before `check_improper()`.", call. = FALSE)
      }
      return(improper(private$.estimate))
    },

#' @description
#' Get a specific estimate from the fitted SemFC model.
#'
#' @param estimate Character string specifying the estimate component to select. 
#'   Possible values are: 
#'   
#'   \code{"lambda"}, \code{"std_lambda"}, \code{"omega"}, \code{"std_omega"}, 
#'   \code{"residual_variance"}, \code{"p_tilde"}, 
#'   \code{"beta"}, \code{"gamma"}, \code{"p_exo"}, \code{"p_endo"}, \code{"psi"},
#'   \code{"p_implied"}, \code{"sigma_implied"}, \code{"r2"}, \code{"T_LS"}, 
#'   \code{"theta"}, and \code{"effect"}. 
#'   
#'   Use \code{"all"} to retrieve all estimates. 
#'
#' @return The requested estimate. Returns \code{NULL} with a warning if:
#'   \itemize{
#'     \item The model has not been fitted yet
#'     \item The requested estimate is not available in the fitted model
#'   }
#'   When \code{estimate = "all"}, returns a list containing all available 
#'   estimates.
#'
#' @examples
#' data(ECSI)
#'
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
#' # Initialize SemFC object with data, model specification 
#' # and svd-based estimation method.
#'
#' sem_model <- SemFC$new(data = A,
#'                        relation_matrix = C,
#'                        mode = rep("reflective", 5),
#'                        scale = FALSE,
#'                        estimator = "svd")
#'                        
#' # Fit the full model with no inference and goodness-of-fit.                        
#'
#' sem_model$fit()
#' 
#' # Get specific estimates
#' lambda <- sem_model$get_estimate("lambda")
#' beta <- sem_model$get_estimate("beta")
#' r2 <- sem_model$get_estimate("r2")
#'
#' # Get all estimates
#' all_estimates <- sem_model$get_estimate("all")
#'
    get_estimate = function(estimate){
      if (is.null(private$.estimate) || length(private$.estimate) == 0L) {
        warning("Model has not been fitted yet. Returning NULL.", call. = FALSE)
        return(NULL)
      }

      if (estimate == "all") {
        return(private$.estimate)
      } else if (!estimate %in% names(private$.estimate)) {
        warning("Estimate '", estimate, 
                "' not found in fitted model. Returning NULL.", 
                call. = FALSE)
        return(NULL)
      }

      return(private$.estimate[[estimate]])
    },

get_model = function(attribute){
  return(private$.model[[attribute]])
},

get_data = function(attribute){
  return(private$.data[[attribute]])}


  ),

  private = list(
#
# ' @field estimator Character string specifying estimation method ("svd" or "ml")
# ' @field data List containing data-related components:
# '   \itemize{
# '     \item \code{data}: The input data (list of blocks)
# '     \item \code{n_row}: Number of observations
# '     \item \code{cov_S}: Covariance matrix of observed variables
# '     \item \code{S_diag_composites}: List of covariance matrices for formative blocks
# '   }
# ' @field model List containing model specification and parameters:
# '   \itemize{
# '     \item \code{relation_matrix}: Square matrix defining structural relationships
# '     \item \code{mode}: Character vector specifying measurement model types
# '     \item \code{n_blocks}: Number of measurement blocks
# '     \item \code{varnames}: List of variable names for each block
# '     \item \code{block_sizes}: Vector of sizes for each measurement block
# '     \item \code{dag}: Logical indicating if structural model is recursive
# '     \item \code{which_exo_endo}: List identifying exogenous/endogenous variables
# '     \item \code{lengths_theta}: Vector of parameter counts
# '     \item \code{p}: Total number of observed variables
# '     \item \code{q}: Total number of free parameters
# '     \item \code{r}: Number of formative blocks
# '     \item \code{dof}: Degrees of freedom
# '     \item \code{scale}: Logical indicating whether to standardize data
# '     \item \code{bias}: Logical indicating bias correction in covariance estimation
# '   }
# ' @field estimate List containing all estimated model parameters:
# '   \itemize{
# '     \item \code{lambda}: List of loading vectors for each block
# '     \item \code{omega}: List of composite weight vectors for formative blocks
# '     \item \code{beta}: Matrix of structural paths between endogenous variables
# '     \item \code{gamma}: Matrix of structural paths from exogenous to endogenous variables
# '     \item \code{residual_variance}: List of residual variances for each observed variable
# '     \item \code{P_EXO}: Correlation matrix of exogenous latent variables
# '     \item \code{P_ENDO}: Correlation matrix of endogenous latent variables
# '     \item \code{P_IMPLIED}: Implied correlation matrix of all latent variables
# '     \item \code{SIGMA_IMPLIED}: Implied covariance matrix of observed variables
# '     \item \code{std_lambda}: List of standardized loadings for each block
# '     \item \code{std_omega}: List of standardized composite weights for formative blocks
# '     \item \code{psi}: Residual covariance matrix of latent variables
# '     \item \code{R2}: Named vector of R-squared values for endogenous latent variables
# '     \item \code{T_LS}: Least squares fit function value
# '     \item \code{theta}: Numeric vector of all free parameters
# '     \item \code{effect}: List of total and indirect effects between latent variables
# '     \item \code{Ptilde}: First-step correlation matrix estimate (SVD only)
# '   }
# ' @field infer_estimate List containing inference results for all parameters:
# '   \itemize{
# '     \item \code{lambda}: Data frame of loadings with SE, z-values, p-values and CI
# '     \item \code{omega}: Data frame of composite weights with SE, z-values, p-values and CI
# '     \item \code{beta}: Data frame of structural paths (endo) with SE, z-values, p-values and CI
# '     \item \code{gamma}: Data frame of structural paths (exo) with SE, z-values, p-values and CI
# '     \item \code{residual_variance}: Data frame of residual variances with SE, z-values, p-values and CI
# '     \item \code{total_effects}: Data frame of total effects with SE, z-values, p-values and CI
# '     \item \code{indirect_effects}: Data frame of indirect effects with SE, z-values, p-values and CI
# '     \item \code{VCOV}: Variance-covariance matrix of parameter estimates (ML only)
# '     \item \code{vcov_effect}: Variance-covariance matrix of effect estimates (ML only)
# '   }
# ' @field boot_rep Integer number of bootstrap replications
# ' @field gof List containing goodness-of-fit statistics:
# '   \itemize{
# '     \item \code{F}: Value of the fit function at the optimal solution
# '     \item \code{reliability}: Named vector of reliability coefficients (Dillon-Goldstein rho) for reflective blocks
# '     \item \code{chi2}: Chi-square test statistic
# '     \item \code{df}: Degrees of freedom for chi-square test
# '     \item \code{pvalue}: P-value of the chi-square test
# '     \item \code{CFI}: Comparative Fit Index
# '     \item \code{TLI}: Tucker-Lewis Index
# '     \item \code{RMSEA}: Root Mean Square Error of Approximation
# '     \item \code{RMSEA_CI}: 90\% confidence interval for RMSEA
# '     \item \code{SRMR}: Standardized Root Mean Square Residual
# '     \item \code{AIC}: Akaike Information Criterion
# '     \item \code{BIC}: Bayesian Information Criterion
# '     \item \code{SABIC}: Sample-size Adjusted BIC
# '     \item \code{bollen_stine}: Bollen-Stine bootstrap p-value (SVD only)
# '   }


    .estimator      = NULL,
    .data           = list(),
    .model          = list(),
    .estimate       = list(),
    .infer_estimate = NULL,
    .boot_rep       = NULL,
    .gof            = list(),



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
      private$.estimator <- 'svd'
      svd_result <- svdSEM(private$.data$data,
                           private$.model$relation_matrix,
                           private$.model$scale,
                           private$.model$mode,
                           private$.model$bias)

      private$.estimate <- svd_result
      theta_svd <- parameters_svd(lambda = svd_result$lambda,
                                  P_EXO = svd_result$p_exo,
                                  G = svd_result$gamma,
                                  B = svd_result$beta,
                                  P_ENDO = svd_result$p_endo,
                                  residual_variance = svd_result$residual_variance,
                                  S_composites = private$.data$S_diag_composites,
                                  model = private$.model)

      private$.estimate$theta <- theta_svd
      private$.estimate$effect <- compute_effect(private$.estimate$beta, private$.estimate$gamma)
      private$.gof$F <- F1(theta_svd, private$.data$cov_S, private$.model)
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

      boot_out <- svdsem_infer(private$.estimate, B, verbose = verbose)
      private$.infer_estimate <- boot_out$result$infer
      private$.gof$bollen_stine <- boot_out$gof


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

      len_theta <- sum(private$.model$lengths_theta)

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
            private$.estimator <- 'ml'
            private$.estimate$theta
          },
          random = runif(len_theta)
        )
      }

      ml_sol <- mlSEM(initial_params, private$.data$cov_S, private$.model, tol)
      theta_ml <- ml_sol$pars
      private$.estimate <- lvm_ml(x = theta_ml, model = private$.model, jac = F)
      private$.estimate$T_LS <- d_LS(private$.data$cov_S, private$.estimate$sigma_implied)

      var_MVs <- lapply(private$.data$data, function(x) diag(cov2(x, bias = private$.model$bias)))
      std_lambda <- mapply("/", private$.estimate$lambda, lapply(var_MVs, sqrt),  SIMPLIFY = FALSE)
      private$.estimate$std_lambda <- std_lambda
      std_omega <- mapply(function(Sjj, lambda_j) solve(Sjj) %*% lambda_j,
                         private$.data$S_diag_composites, std_lambda[private$.model$mode == "formative"],
                         SIMPLIFY = FALSE)
      names(std_omega) <- names(std_lambda[private$.model$mode == "formative"])
      private$.estimate$std_omega <- std_omega
      private$.estimate$theta <- theta_ml
      private$.estimate$effect <- compute_effect(private$.estimate$beta, private$.estimate$gamma)
      private$.gof$F <- F1(theta_ml, private$.data$cov_S, private$.model)


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
      theta_ml <- private$.estimate$theta
      S <- private$.data$cov_S
      N <- private$.data$n_row

      ml_infer_estimate <- mlSEM_infer(theta_ml, S, private$.model, N, private$.estimate)


      private$.infer_estimate <- ml_infer_estimate$estimate
      private$.infer_estimate$VCOV <- ml_infer_estimate$VCOV
      private$.infer_estimate$vcov_effect <- ml_infer_estimate$vcov_effect
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

      estimator <- private$.estimator

      res_gof <- list()

      # reliability only for relflective block (Dillon)
      if (sum(private$.model$mode == "reflective") > 0){
        res_reliability <- reliability('Dillon', private$.estimate$lambda, private$.estimate$residual_variance)
        res_gof$reliability <- res_reliability
      }


      if (estimator == 'svd' && is.null(private$.gof$bollen_stine)){
        bollen_stine <- svdSEM_gof(private$.estimate, B)
        res_gof$bollen_stine <- bollen_stine
      } else if (estimator == 'ml'){
        p <- private$.model$p
        q <- private$.model$q
        r <- private$.model$r
        F <- private$.gof$F
        N <- private$.data$n_row
        S <- private$.data$cov_S
        Sigma <- private$.estimate$sigma_implied

        res_gof <- c(res_gof, semML_gof(p, q, r, F, N, S, Sigma))

      }
      private$.gof <- c(private$.gof, res_gof)

    }
  )
)


