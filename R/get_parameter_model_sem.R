
#' Extract Parameters for Structural Equation Models
#'
#' This function processes input data and computes various parameters required
#' for structural equation modeling (SEM), including covariance matrices, block sizes,
#' and variable names. It also identifies composite covariance blocks for formative models.
#'
#' @param data A list of data frames or matrices, where each element represents
#'   a block of observed variables.
#' @param mode A character vector specifying the measurement model type for each block:
#'   "formative" or "reflective".
#'
#' @return A list containing:
#'   \item{data}{The input data (list of blocks).}
#'   \item{n_blocks}{The number of blocks in the data.}
#'   \item{n_row}{The number of rows (observations) in the combined data.}
#'   \item{varnames}{A list of variable names for each block.}
#'   \item{block_sizes}{A vector containing the number of variables in each block.}
#'   \item{S}{The covariance matrix of the combined data.}
#'   \item{S_diag_composites}{A list of covariance matrices for formative blocks.}
#'
#' @details The function calculates the covariance matrix of the combined data
#'   and extracts diagonal blocks corresponding to each set of observed variables.
#'   For blocks specified as "formative", their covariance matrices are stored separately.
#'
#' @examples
#' \dontrun{
#' data <- list(
#'   block1 = matrix(rnorm(100), ncol = 5),
#'   block2 = matrix(rnorm(80), ncol = 4)
#' )
#' mode <- c("formative", "reflective")
#' parameters <- get_parameter_model_sem(data, mode)
#' print(parameters$S_diag_composites)
#' }
#'
#'
#'
#' @importFrom igraph graph_from_adjacency_matrix is_dag
#'
#'
#'
get_parameter_model_sem <- function(data, mode, relation_matrix){

  X <- do.call(cbind, data)
  S <- cov(X)
  n_row <- NROW(X)

  block_sizes <- sapply(data, NCOL)
  n_blocks <- length(block_sizes)
  varnames <- lapply(data, function(x) colnames(x))


  # get composite covariance bloc

  start_indices <- unname(cumsum(c(1, head(block_sizes, -1))))
  end_indices <- unname(cumsum(block_sizes))
  S_diag <- mapply(function(start, end) {
    S[start:end, start:end]
  }, start_indices, end_indices, SIMPLIFY = FALSE)

  S_diag_composites <- S_diag[mode == 'formative']

  dag <- igraph::is_dag(graph_from_adjacency_matrix(relation_matrix))
  which_exo_endo <- ind_exo_endo(relation_matrix)
  lengths_theta <- get_lengths_theta(which_exo_endo, block_sizes, mode, dag)

  p <- sum(block_sizes)
  q <- sum(lengths_theta)
  r <- sum(mode == "formative")

  dof  <- (p * (p+1)/2) - q + r

  out <- list(
    data = data,
    n_blocks = n_blocks,
    n_row = n_row,
    varnames = varnames,
    block_sizes = block_sizes,
    cov_S = S,
    S_diag_composites = S_diag_composites,
    dag = dag,
    which_exo_endo = which_exo_endo,
    lengths_theta = lengths_theta,
    p = p,
    q = q,
    r = r,
    dof = dof

  )


  return(out)




}