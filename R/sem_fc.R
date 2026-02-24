
#' @import R6
#' @name SemFC
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
    #' Fit the complete model with inference and goodness-of-fit
    #'
    #' @param infer Logical indicating whether to perform statistical inference (default: FALSE)
    #' @param B Integer number of bootstrap replications for svd (default: 1000)
    #' @param initialization Character string or numeric vector specifying initialization method
    #'   for ML estimation. Ignored when estimator is "svd".
    #'   \itemize{
    #'     \item \code{"svd"} (default): Use SVD estimates as starting values
    #'     \item \code{"random"}: Use random starting values
    #'     \item \code{numeric vector}: Use provided values as starting parameters
    #'       (must have length equal to total number of model parameters)
    #'   }
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
    #' model$fit(estimator = "ml", B = 500, initialization = "svd")
    #' }

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
    #' Extract parameter estimates from the fitted model
    #'
    #' @param standardized Logical indicating whether to include standardized
    #'   estimates (default: FALSE). When TRUE, adds a `std.all` column with
    #'   fully standardized coefficients.
    #'
    #' @details
    #' Returns a data frame containing all estimated parameters including:
    #' \itemize{
    #'   \item \code{lambda}: Loadings (measurement model)
    #'   \item \code{omega}: Composite weights (for formative blocks)
    #'   \item \code{beta}: Structural paths between endogenous variables
    #'   \item \code{gamma}: Structural paths from exogenous to endogenous variables
    #'   \item \code{residualvariance}: Residual variances of observed variables
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
    #' \dontrun{
    #' model$fit(infer = TRUE)
    #' model$parameterEstimates()
    #' model$parameterEstimates(standardized = TRUE)
    #' }

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
    #' Check for improper solutions in the estimated model
    #'
    #' @details
    #' Detects types of inadmissible or improper solutions including:
    #' \itemize{
    #'   \item \code{reliability_coef}: Reliability coefficients for each block
    #'   \item \code{P_IMPLIED}: Implied correlation matrix of latent variables
    #'   \item \code{residual_variance}: List of residual variances
    #'   \item \code{std_lambda}: Standardized loadings
    #'   \item \code{SIGMA_IMPLIED}: Implied covariance matrix of observed variables
    #'   \item \code{R2}: R-squared values for endogenous latent variables
    #'   \item \code{psi}: Residual covariance matrix of latent variables
    #'   \item \code{Ptilde}: First estimation of correlation matrix of latent variables (only for SVDSEM)
    #' }
    #'
    #'
    #' @return Named logical vector of length 8 (or 9 for SVDSEM) indicating presence of each type
    #'   of improper solution:
    #' \describe{
    #'   \item{RELIABILITY_COEF}{\code{TRUE} if any reliability coefficient is outside (0, 1).}
    #'   \item{RHO_JH}{\code{TRUE} if any correlation in P_IMPLIED is outside (-1, 1).}
    #'   \item{P_IMPLIED}{\code{TRUE} if P_IMPLIED has negative eigenvalues (not positive definite).}
    #'   \item{THETA_JH}{\code{TRUE} if any residual variance is negative.}
    #'   \item{STD_LAMBDA}{\code{TRUE} if any standardized loading is outside (-1, 1).}
    #'   \item{SIGMA_IMPLIED}{\code{TRUE} if SIGMA_IMPLIED has negative eigenvalues (not positive definite).}
    #'   \item{R2}{\code{TRUE} if any R-squared is outside (0, 1).}
    #'   \item{PSI}{\code{TRUE} if PSI has negative eigenvalues (not positive definite).}
    #'   \item{P_TILDE}{\code{TRUE} if P_TILDE has negative eigenvalues (only for SVDSEM).}
    #' }
    #'
    check_improper = function(){
      return(improper(self$estimate))
    }
  ),

  private = list(

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
    #
    # ' @description
    # ' Perform statistical inference for SVD estimates using bootstrap
    # '
    # ' @param B Integer number of bootstrap replications (default: 1000)
    # ' @param verbose Logical indicating whether to print progress messages (default: TRUE)
    # '
    # ' @details
    # ' Uses non-parametric bootstrap to estimate standard errors, confidence intervals,
    # ' and p-values for all model parameters.
    # '
    # ' @return Invisible self (for method chaining)
    # ' @keywords internal
    svd_infer = function(B = 1000, verbose = TRUE){

      boot_out <- svdsem_infer(self$estimate, B, verbose = TRUE)
      self$infer_estimate <- boot_out$result$infer
      self$gof$bollen_stine <- boot_out$gof


    },

    # ' @description
    # ' Fit the model using Maximum Likelihood (ML) estimation
    # '
    # ' @param initialization Character string or numeric vector specifying initialization method:
    # '   \itemize{
    # '     \item \code{"svd"} (default): Use SVD estimates as starting values
    # '     \item \code{"random"}: Use random starting values
    # '     \item \code{numeric vector}: Use provided values as starting parameters
    # '       (must have length equal to total number of model parameters)
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

      var_MVs <- lapply(self$data$data, function(x) diag(cov2(x, bias = self$model$bias)))
      std_lambda <- mapply("/", self$estimate$lambda, lapply(var_MVs, sqrt),  SIMPLIFY = FALSE)
      self$estimate$std_lambda <- std_lambda

      self$estimate$theta <- theta_ml
      self$estimate$effect <- compute_effect(self$estimate$beta, self$estimate$gamma)
      self$gof$F <- F1(theta_ml, self$data$cov_S, self$model)


    },


    # ' @description
    # ' Perform asymptotic statistical inference for ML estimates
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
    # ' @param B Integer number of bootstrap replications for Bollen-Stine test (default: 1000).
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

        res_gof <- c(res_gof, semML_gof(p, q, r, F, N, S, Sigma))

      }
      self$gof <- c(self$gof, res_gof)

    }




  )
)


