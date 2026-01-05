
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
#' @param relation_matrix A square adjacency matrix representing the relationships
#'  among latent variables in the SEM.
#' @param bias Logical indicating whether to use biased covariance estimation.
#'
#' @return A list with two main components:
#' \describe{
#'   \item{data}{A list containing:
#'     \itemize{
#'       \item \code{data}: The input data (list of blocks).
#'       \item \code{n_row}: The number of rows (observations) in the combined data.
#'       \item \code{cov_S}: The covariance matrix of the combined data.
#'       \item \code{S_diag_composites}: A list of covariance matrices for formative blocks.
#'     }
#'   }
#'   \item{model}{A list containing:
#'     \itemize{
#'       \item \code{relation_matrix}: The input relation matrix.
#'       \item \code{mode}: The input mode vector.
#'       \item \code{n_blocks}: The number of blocks in the data.
#'       \item \code{varnames}: A list of variable names for each block.
#'       \item \code{block_sizes}: A vector containing the number of variables in each block.
#'       \item \code{dag}: Logical indicating if the structural model is a DAG.
#'       \item \code{which_exo_endo}: List containing indices and structure information for
#'       exogenous and endogenous latent variables (output from `ind_exo_endo()`).
#'       \item \code{lengths_theta}: Vector of estimation parameter counts.
#'       \item \code{lengths_cov_parameter}: Vector of covariance estimation parameter counts for each block.
#'       \item \code{p}: Number of observed variables in the model
#'       \item \code{q}: Number of estimated parameters in the model
#'       \item \code{r}: Number of formative blocks.
#'       \item \code{dof}: Degrees of freedom for model testing.
#'     }
#'   }
#' }
#'
#' @details The function calculates the covariance matrix of the combined data
#'   and extracts diagonal blocks corresponding to each set of observed variables.
#'   For blocks specified as "formative", their covariance matrices are stored separately.
#'
#'
#'
#'
#' @importFrom igraph graph_from_adjacency_matrix is_dag
#'
#'
#' @keywords internal
get_parameter_model_sem <- function(data, mode, relation_matrix, bias){

  X <- do.call(cbind, data)
  S <- cov2(X, bias = bias)
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
  # lengths of covariance parameters for each block
  lengths_cov_parameter <- ifelse(mode == "formative", (block_sizes^2 + block_sizes) / 2, block_sizes)

  p <- sum(block_sizes)
  q <- sum(lengths_theta)
  r <- sum(mode == "formative")

  dof  <- (p * (p+1)/2) - q + r

  out <- list(
    data = list(
      data = data,
      n_row = n_row,
      cov_S = S,
      S_diag_composites = S_diag_composites
    ),


    model = list(
      relation_matrix = relation_matrix,
      mode = mode,
      n_blocks = n_blocks,
      varnames = varnames,
      block_sizes = block_sizes,
      dag = dag,
      which_exo_endo = which_exo_endo,
      lengths_theta = lengths_theta,
      lengths_cov_parameter = lengths_cov_parameter,
      p = p,
      q = q,
      r = r,
      dof = dof
    )

  )


  return(out)




}