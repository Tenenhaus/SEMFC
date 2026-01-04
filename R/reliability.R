#' Compute Reliability for Formative Blocks
#'
#' This function calculates the reliability of formative measurement blocks
#' using the specified metric. Currently, only the "Dillon" metric is supported.
#'
#' @param metric A character string specifying the reliability metric to use.
#'   Currently, only "Dillon" is implemented.
#' @param lambdas A list of factor loadings for each block.
#' @param residual_variances A list of residual variances for each block.
#'
#' @return A named numeric vector containing the reliability values for each block.
#'
#' @details The Dillon reliability metric is calculated as the squared sum of
#'   loadings divided by the sum of the squared loadings and the residual variances.
#'
#' @examples
#' \dontrun{
#' lambdas <- list(block1 = c(0.8, 0.7), block2 = c(0.9, 0.6, 0.5))
#' residuals <- list(block1 = c(0.2, 0.3), block2 = c(0.1, 0.4, 0.5))
#' reliability_values <- reliability(
#'  metric = "Dillon",
#'  lambdas = lambdas,
#'  residual_variances = residuals
#' )
#' print(reliability_values)
#' }
#' @export

# pour bloc formatif
reliability <- function (metric, lambdas, residual_variances){
  if (metric=='Dillon'){
    # compute sum of loadings and residual variances for each block
    sum_loadings <- sapply(lambdas, function(x) sum(as.numeric(x)))
    sum_residuals <- sapply(residual_variances, function(x) sum(as.numeric(x)))

    #remove formative blocks (blocks without residual variances)
    sum_loadings <- sum_loadings[names(sum_loadings) %in% names(sum_residuals)]

    # compute numerator and denominator for Dillon's rho
    numerator <- sum_loadings^2
    denominator <- numerator + sum_residuals

    # compute reliability results
    results <- numerator / denominator

    return(results)
  }
}

