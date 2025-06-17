

get_parameter_model_sem <- function(data, mode){

  X <- do.call(cbind, data)
  S <- cov(X)
  n_row <- NROW(X)

  block_sizes <- sapply(data, NCOL)
  n_blocks <- length(block_sizes)


  # get composite covariance bloc

  start_indices <- unname(cumsum(c(1, head(block_sizes, -1))))
  end_indices <- unname(cumsum(block_sizes))
  S_diag <- mapply(function(start, end) {
    S[start:end, start:end]
  }, start_indices, end_indices, SIMPLIFY = FALSE)

  S_diag_composites <- S_diag[mode == 'formative']

  out <- list(
    data = data,
    n_blocks = n_blocks,
    n_row = n_row,
    block_sizes = block_sizes,
    S = S,
    S_diag_composites = S_diag_composites

  )


  return(out)












}