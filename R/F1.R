
#' Compute Log-Likelihood for ML Estimation of Structural Equation Model
#'
#' This function computes the log-likelihood objective function for maximum
#' likelihood estimation of a structural equation model with latent variables.
#' It evaluates the discrepancy between the sample covariance matrix and the
#' model-implied covariance matrix.
#'
#' @param x Numeric vector containing all model parameters (loadings, correlations,
#'   path coefficients, and variance/covariance parameters).
#' @param S Sample covariance matrix of observed variables.
#' @param model A list containing model specifications with the following elements:
#'   \describe{
#'     \item{block_sizes}{Integer vector specifying the number of indicators in each block.}
#'     \item{mode}{Character vector indicating the measurement mode for each block
#'       ("formative" or "reflective").}
#'     \item{lengths_parameter}{Integer vector specifying the length of each parameter
#'       group in `x` (loadings, exogenous correlations, gamma, beta, endogenous correlations,
#'       covariance of composite blocks and residual variances).}
#'     \item{which_exo_endo}{List containing indices and structure information for
#'       exogenous and endogenous latent variables (output from `ind_exo_endo()`).}
#'    \item{dag}{Logical indicating whether the structural model is recursive (FALSE)
#'       or non-recursive (TRUE).}
#'   }
#'
#' @return Numeric scalar representing the log-likelihood value to be minimized.
#'
#' @details
#' The log-likelihood function is computed as:
#' \deqn{\log|\Sigma(\theta)| + \text{tr}(S\Sigma(\theta)^{-1}) - \log|S| - p}
#' where:
#' - \eqn{\Sigma(\theta)} is the model-implied covariance matrix
#' - \eqn{S} is the sample covariance matrix
#' - \eqn{p} is the number of observed variables
#'
#' The function first computes the implied covariance matrix using `lvm_ml()`,
#' then evaluates the fit between sample and implied covariances.
#'
#' @seealso \code{\link{lvm_ml}} for implied covariance matrix computation,
#'   \code{\link{ind_exo_endo}} for structure information generation
#'
#'
#' @examples
#' # Example with 6 blocks: 4 formative exogenous, 2 reflective endogenous
#' set.seed(123)
#' n <- 300
#'
#'
#' # Generate sample data (6 blocks with 3 indicators each)
#' X <- matrix(rnorm(n * 18), ncol = 18)
#' S <- cov(X)
#'
#' # Model specification
#' block_sizes <- c(3, 3, 3, 3, 3, 3)
#' mode <- c("formative", "formative", "formative", "formative",
#'           "reflective", "reflective")
#'
#' # Structure information (hard-coded for this example)
#' which_exo_endo <- list(
#'   Hi = list(
#'     c(LV1 = 1, LV2 = 2),
#'     c(LV3 = 3, LV4 = 4)
#'   ),
#'   Ji = list(
#'     c(LV6 = 6),
#'     c(LV5 = 5)
#'   ),
#'   ind_endo = c(LV5 = 5, LV6 = 6),
#'   ind_exo = c(LV1 = 1, LV2 = 2, LV3 = 3, LV4 = 4)
#' )
#'
#' dag <- FALSE
#'
#' varnames <- list(
#'  LV1 = c("X11", "X12", "X13"),
#'  LV2 = c("X21", "X22", "X23"),
#'  LV3 = c("X31", "X32", "X33"),
#'  LV4 = c("X41", "X42", "X43"),
#'  LV5 = c("X51", "X52", "X53"),
#'  LV6 = c("X61", "X62", "X63")
#')
#'
#'
#' # Parameter vector structure:
#' # - 18 loadings (3 per  block)
#' # - 6 exogenous correlations
#' # - 4 non zero gamma coefficients
#' # - 2 non zero beta coefficient
#' # - 1 endogenous correlations
#' # - 30 covariances and residual variances
#' lengths_theta <- c(18, 6, 4, 2, 1, 30)
#'
#'
#' model <- list(
#'  block_sizes = block_sizes,
#' mode = mode,
#' varnames = varnames,
#' lengths_theta = lengths_theta,
#' which_exo_endo = which_exo_endo,
#' dag = dag
#' )
#'
#'
#'
#' # Initialize parameters
#' x  <- rnorm(61)
#' # Compute log-likelihood
#' f <- F1(x, S, model)
#' print(f)
#'
#' @export
F1 <- function(x, S, model){


  implied_S <- lvm_ml(x, model, jac = FALSE)$SIGMA_IMPLIED

  # check for positive definiteness

  eigvals <- eigen(implied_S, symmetric = TRUE, only.values = TRUE)$values
  if (any(eigvals <= 0.0001)){
    return(1e+5)
  }

  opt <- log(det(implied_S)) +
    sum(diag(S%*%solve(implied_S))) -
    log(det(S)) -
    NCOL(S)

  return(opt)

}