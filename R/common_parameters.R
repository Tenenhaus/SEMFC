#' Construct Block Diagonal Covariance Structures for Formative and Reflective Blocks
#'
#' This function creates a list of block diagonal covariance matrix for both formative and reflective blocks
#' using the parameter vector `x`, block sizes, and the starting index of covariance parameters.
#'
#' @param x A numeric vector of parameters from which covariance or variance values are extracted.
#' @param mode A character vector indicating the mode ("formative" or "reflective") for each block.
#' @param block_sizes A numeric vector where each element specifies the number of indicators
#'   of each block.
#' @param initial_start_index_cov An integer indicating the starting index in `x` for covariance/variance
#'   parameters.
#'
#' @return A list of block matrices, where each block corresponds to either:
#'   - A covariance matrix for formative blocks.
#'   - A diagonal matrix for reflective blocks.
#'
#' @examples
#' \dontrun{
#' set.seed(27)
#' x <- rnorm(61)
#' mode <- c("formative", "formative", "formative", "formative", "reflective", "reflective")
#' block_sizes <- rep(3, 6)
#' initial_start_index_cov <- 32
#' get_bdiag(x, mode, block_sizes, initial_start_index_cov)}
#' @keywords internal

get_bdiag <- function(x, mode, block_sizes, initial_start_index_cov) {


  # number of blocks
  J <- length(block_sizes)
  # list of lengths of the upper values in the cov matrix for the composite block i or
  # of the diagonal values for formative for each block
  lengths_values_cov <- block_sizes
  lengths_values_cov[mode == "formative"] <- (block_sizes[mode == "formative"]^2 + block_sizes[mode == "formative"]) / 2
  # number of parameters for covariance
  total_cov_parameter <- sum(lengths_values_cov)
  end_endex_cov <- initial_start_index_cov + total_cov_parameter - 1

  # part of the vector corresponding to covariance blocks
  extracted_parameters_cov <- x[initial_start_index_cov:end_endex_cov]
  # list of parameters corresponding to each covariance bloc
  list_cov <- split(extracted_parameters_cov,
                    rep(seq_along(lengths_values_cov), lengths_values_cov))
  names(list_cov) <- names(lengths_values_cov)

  # Building list of formative matrices
  S_composites <- lapply(list_cov[mode == "formative"], build_formative_S_diag)

  # Building list of reflectives matrices
  reflective_blocks <- lapply(list_cov[mode == "reflective"], diag)

  BDIAG <- vector("list", J)
  BDIAG[mode == "formative"] <- S_composites
  BDIAG[mode == "reflective"] <- reflective_blocks
  names(BDIAG) <- names(lengths_values_cov)

  return(BDIAG)

}








#' Construct Correlation Matrix from Parameter Vector
#'
#' This function constructs a symmetric correlation matrix for a set of latent variables
#' based on the parameter vector `x` and the starting index for correlation-related parameters.
#'
#' @param x A numeric vector of parameters. Correlation-related values are extracted starting
#'   from the index specified by `start_index`.
#' @param latent_variables A vector of names or identifiers for the latent variables, used
#'   to determine the dimensions of the correlation matrix.
#' @param start_index An integer specifying the index in `x` where the correlation parameters
#'   begin.
#'
#' @return A symmetric correlation matrix (dim x dim), where `dim` is the number of latent variables.
#'
#' @examples
#' \dontrun{
#' latent_variables <- c("LV1", "LV2", "LV3", "LV4")
#' x <- c(0.5, 0.3, 0.2, 0.5, 0.4, 0.2) # lower triangle values for a 4x4 matrix
#' start_index <- 1
#' get_correlation_coeff(x, latent_variables, start_index)}
#' @keywords internal
#'
get_correlation_coeff <- function(x, latent_variables, start_index) {


  # dim  matrix (number of exogeneous or endogeneous latent variables)
  dim <- length(latent_variables)
  #length lower values vector from this symetric matrix:
  length_correl <- dim * (dim - 1) / 2
  end_index <- start_index + length_correl - 1
  lower_values <- x[start_index:end_index]

  # Building of the correlation matrix
  P_matrix <- diag(1, dim, dim)
  if (dim > 1) {
    P_matrix[lower.tri(P_matrix)] <- lower_values
    P_matrix[upper.tri(P_matrix)] <- t(P_matrix)[upper.tri(P_matrix)]
  }

  return(P_matrix)
}


#' Extract Loadings from Parameter Vector
#'
#' This function retrieves the loadings for each block of observed variables
#' based on the parameter vector `x` and the input block sizes.
#'
#' @param x Numeric vector containing all model parameters (loadings, correlations,
#'   path coefficients, and variance/covariance parameters).
#' @param block_sizes A numeric vector where each element represents the number
#'   of observed variables in the corresponding block of latent variables.
#'
#' @return A list of numeric vectors, where each vector contains the loadings
#'   for the corresponding block of manifest variables.
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' block_sizes <- c(3, 3, 3, 3, 3, 3)
#' x  <- rnorm(61)
#' get_loadings(x, block_sizes)}
#' @keywords internal
get_loadings <- function(x, block_sizes) {


  # total number of emergent variables
  y <- sum(block_sizes)
  loadings <- setNames(
    lapply(seq_along(block_sizes), function(b) {
      setNames(
        split(x[1:y], rep(seq_along(block_sizes), block_sizes))[[b]],
        paste0(names(block_sizes[b]), 1:block_sizes[b])
      )
    }),
    names(block_sizes)
  )


  return(loadings)

}


#' Construct Path Coefficient Matrix from Parameter Vector
#'
#' This function generates a path coefficient matrix for endogenous and exogenous variables
#' based on the parameter vector `x` and the provided relationships between variables.
#'
#' @param x A numeric vector of parameters. Contains path coefficient values to be extracted.
#' @param row_exo_endo A named numeric vector representing the row indices (endogenous variables).
#' @param col_exo_endo A named numeric vector representing the column indices (exogenous or endogenous variables).
#' @param C A matrix defining the relationships between variables (1 indicates a relationship exists).
#' @param initial_start_index An integer specifying the index in `x` where the path coefficient
#'   parameters start.
#'
#' @return A matrix where rows correspond to endogenous variables and columns to exogenous or
#'   endogenous variables, with path coefficients filled in where applicable.
#'
#' @examples
#' \dontrun{
#' x <- c(0.8, 0.5, 0.3, 0.7, 0.8, 0.4)
#' row_exo_endo <- c(X5 = 5, X6 = 6)
#' col_exo_endo <- c(X1 = 1, X2 = 2, X3 = 3, X4 = 4)
#' initial_start_index <- 1
#' C <- matrix(c(0, 0, 0, 0, 1, 0,
#'             0, 0, 0, 0, 1, 0,
#'             0, 0, 0, 0, 0, 1,
#'             0, 0, 0, 0, 0, 1,
#'             0, 0, 0, 0, 0, 1,
#'             0, 0, 0, 0, 1, 0), 6, 6, byrow = TRUE)
#' get_path_coeff(x, row_exo_endo, col_exo_endo, C, initial_start_index)}
#' @keywords internal

get_path_coeff <- function(x, row_exo_endo, col_exo_endo, C, initial_start_index) {


  vec_m <- as.vector(t(C[col_exo_endo, row_exo_endo, drop = FALSE]))
  len_param <- sum(vec_m)
  if (len_param == 0){
    return(matrix(0, nrow = length(row_exo_endo), ncol = length(col_exo_endo)))
  }
  end_index <- initial_start_index + len_param - 1
  M_select <- diag(length(vec_m))[, vec_m == 1, drop = FALSE]
  vec_regression <- M_select%*%x[initial_start_index:end_index]

  Matrix_path <- matrix(vec_regression, nrow = length(row_exo_endo), ncol = length(col_exo_endo))

  return(Matrix_path)
}



#' Construct Symmetric Formative S Matrix
#'
#' This function builds a symmetric covariance matrix (S_composite_i) for a formative block,
#' using the lower triangular values provided in `lower_values_composite_i`.
#'
#' @param lower_values_composite_i A numeric vector containing the values for the lower
#'   triangle (including the diagonal) of the symmetric matrix.
#'
#' @return A symmetric formative covariance matrix.
#'
#' @examples
#' \dontrun{
#' upper_values <- c(1, 0.3, 0.5, 1, 0.2, 1)
#' build_formative_S_diag(upper_values)}
#' @keywords internal


build_formative_S_diag <- function(lower_values_composite_i) {


  # dim of the matrix extracted from the length of the lower values
  dim_block <- (-1 + sqrt(1 + 8 * length(lower_values_composite_i))) / 2
  S_composite_i <- matrix(0, dim_block, dim_block)
  S_composite_i[lower.tri(S_composite_i, diag = TRUE)] <- lower_values_composite_i
  S_composite_i[upper.tri(S_composite_i)] <- t(S_composite_i)[upper.tri(S_composite_i)]

  return(S_composite_i)
}






#' Convert SVD-SEM Parameters to Vector Form
#'
#' Transforms structural equation model parameters from matrix/list format into
#' a single concatenated parameter vector for optimization or inference. The function
#' handles both reflective and formative measurement models.
#'
#' @param lambda List of loading vectors, one per block. Each element contains the
#'   loadings for the indicators in that block.
#' @param P_EXO Symmetric correlation/covariance matrix of exogenous latent variables.
#' @param G Matrix of gamma coefficients (exogenous → endogenous paths). Rows represent
#'   endogenous variables, columns represent exogenous variables.
#' @param B Matrix of beta coefficients (endogenous → endogenous paths). Rows and columns
#'   represent endogenous variables.
#' @param P_ENDO Symmetric correlation/covariance matrix of endogenous latent variables.
#' @param residual_variance List of residual variances for reflective blocks. Each element
#'   contains the residual variances for the indicators in that block.
#' @param S_composites List of empirical covariance matrices for formative composite blocks.
#'   Each element is a covariance matrix for a formative block.
#' @param model A list containing the model structure, including:
#'   \describe{
#'     \item{n_blocks}{Number of blocks in the model.}
#'     \item{mode}{Character vector indicating the measurement mode for each block:
#'       \code{"reflective"} or \code{"formative"}.}
#'     \item{dag}{Logical; if \code{TRUE}, endogenous correlations are excluded
#'       (directed acyclic graph assumption).}
#'   }
#'
#' @return Unnamed numeric vector containing all free parameters in the following order:
#'   \enumerate{
#'     \item Loadings (all blocks concatenated)
#'     \item Exogenous correlations (lower triangle of P_EXO, excluding diagonal)
#'     \item Gamma coefficients (non-zero elements, column by column)
#'     \item Beta coefficients (non-zero elements, column by column)
#'     \item Endogenous correlations (lower triangle of P_ENDO, excluding diagonal)
#'     \item Residual variances (reflective blocks) or composite covariances
#'       (lower triangle with diagonal, formative blocks)
#'   }
#'
#' @details
#' For formative blocks, the lower triangular part (including diagonal) of the
#' empirical covariance matrix is included. For reflective blocks, only residual
#' variances are included.
#'
#'
#' @keywords internal

parameters_svd <- function(lambda,
                           P_EXO,
                           G,
                           B,
                           P_ENDO,
                           residual_variance,
                           S_composites,
                           model){

  J <- model$n_blocks
  mode <- model$mode
  dag <- model$dag

  S_composites_lower <- lapply(S_composites, function(mat) mat[lower.tri(mat, diag = TRUE)])

  # empirical covariance for composite blocks or residual_variance for reflective
  diag_jj <- vector("list", J)
  diag_jj[mode == "formative"] <- S_composites_lower
  diag_jj[mode != "formative"] <- residual_variance


  vect_lambda <- Reduce("c", lambda)
  vect_exo <- P_EXO[lower.tri(P_EXO)]
  vect_endo <- P_ENDO[lower.tri(P_ENDO)]
  if (dag){
      vect_endo <- numeric(0)
  }

  vect_beta <- as.vector(B)[B !=0]
  vect_gamma <- as.vector(G)[G !=0]
  vect_cov <- Reduce("c", diag_jj)



  theta_vect <-
    Reduce("c",
      c(vect_lambda,
        vect_exo,
        vect_gamma,
        vect_beta,
        vect_endo,
        vect_cov
      )
    )



  return(unname(theta_vect))
}





compute_svd_initialization <- function(data, model) {

  svd_result <- svdSEM(
    data$data,
    model$relation_matrix,
    model$scale,
    model$mode,
    model$bias
  )

  theta_svd <- parameters_svd(
    lambda = svd_result$lambda,
    P_EXO = svd_result$p_exo,
    G = svd_result$gamma,
    B = svd_result$beta,
    P_ENDO = svd_result$p_endo,
    residual_variance = svd_result$residual_variance,
    S_composites = data$S_diag_composites,
    model = model
  )

  return(theta_svd)
}
