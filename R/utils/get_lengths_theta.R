get_lengths_theta <- function(C, block_sizes, mode){

  which_exo_endo <- ind_exo_endo(C)
  n <- which_exo_endo$ind_exo
  m <- which_exo_endo$ind_endo

  ##########################################################
  ####### number of parameter for each part ################
  ##########################################################
  number_loadings <- sum(block_sizes)
  number_upper_values_exo <- length(n) * (length(n) - 1) / 2
  # number_non_zero_G <- sum(lengths(which_exo_endo$Hi))
  number_non_zero_G <- sum(length(unlist(which_exo_endo$Hi)[unlist(which_exo_endo$Hi) != 0]))
  # number_non_zero_B <- sum(lengths(which_exo_endo$Ji))
  number_non_zero_B <- sum(length(unlist(which_exo_endo$Ji)[unlist(which_exo_endo$Ji) != 0]))
  number_upper_values_endo <- length(m) * (length(m) - 1) / 2
  number_cov <- sum(ifelse(mode == "formative", (block_sizes^2 + block_sizes) / 2, block_sizes))

  lengths_parameter <- c(number_loadings,
                         number_upper_values_exo,
                         number_non_zero_G,
                         number_non_zero_B,
                         number_upper_values_endo,
                         number_cov)

  return(lengths_parameter)

}