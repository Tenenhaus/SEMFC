

# compute the length of theta model parameter

length_parameters_svd <- function (block_sizes, which_exo_endo){
  n <- which_exo_endo$ind_exo
  m <- which_exo_endo$ind_endo

  length_lambda <- sum(block_sizes)
  length_gamma <- sapply(which_exo_endo$Hi,
                       function(sublist) {
                          ifelse((length(sublist) == 1 && sublist[[1]] == 0),
                                 return(0),
                                 return(length(sublist)))
                       })
  length_beta <- sapply(which_exo_endo$Ji,
                     function(sublist) {
                        ifelse((length(sublist) == 1 && sublist[[1]] == 0),
                               return(0),
                               return(length(sublist)))
                     })

  dim_exo <- length(n)
  length_exo <- dim_exo * (dim_exo - 1) / 2

  dim_endo <- length(m)
  length_endo <- dim_endo * (dim_endo - 1) / 2

  length_residual_variance <- sum(block_sizes[mode == 'reflective'])

  length_cov_composite <- sum(unlist(lapply(block_sizes[mode == 'formative'], function(n) n * (n + 1) / 2)))

  len_parameters <- length_lambda +
    length_gamma + length_beta +
    length_exo + length_endo +
    length_residual_variance + length_cov_composite

  return(len_parameters)

}
