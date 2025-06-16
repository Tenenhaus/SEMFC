#########################################
#       return \Sigma(\theta) for       #
#         hessian computation           #
#########################################

s_implied <- function(x, S, jac = TRUE){
  #loadings
  l1 = x[1:3] ; l2 = x[4:6] ; l3 = x[7:9]
  l4 = x[10:12] ; l5 = x[13:15] ;   l6 = x[16:18]
  L = bdiag(list(l1, l2, l3, l4, l5, l6))
  
  # correlations between exogeneous
  P_EXO = matrix(c(1, x[19], x[20], x[22],
                   x[19], 1, x[21], x[23],
                   x[20], x[21], 1, x[24],
                   x[22], x[23], x[24], 1), 4, 4)
  
  # path coefficients
  G = matrix(c(x[25], x[26], 0, 0,
               0, 0, x[27], x[28]), 2, 4, byrow = TRUE)
  
  # path coefficients
  B = matrix(c(0, x[29],
               x[30], 0), 2, 2, byrow = TRUE)
  
  # Correlations between endogeneous
  P_ENDO = matrix(c(1, x[31],
                    x[31], 1), 2, 2)
  
  # Correlations between Latent/emergent variables
  R = rbind(cbind(P_EXO, P_EXO%*%t(G)%*%t(solve(diag(2) - B))),
            cbind(solve(diag(2) - B)%*%G%*%P_EXO, P_ENDO))
  
  #cov between MVs
  S1 = matrix(c(x[32], x[33], x[35],
                x[33], x[34], x[36],
                x[35], x[36], x[37]), 3, 3)
  #cov between MVs
  S2 = matrix(c(x[38], x[39], x[41],
                x[39], x[40], x[42],
                x[41], x[42], x[43]), 3, 3)
  #cov between MVs
  S3 = matrix(c(x[44], x[45], x[47],
                x[45], x[46], x[48],
                x[47], x[48], x[49]), 3, 3)
  #cov between MVs
  S4 = matrix(c(x[50], x[51], x[53],
                x[51], x[52], x[54],
                x[53], x[54], x[55]), 3, 3)
  
  T5 = diag(x[56:58])
  T6 = diag(x[59:61])
  
  implied_S = L%*%R%*%t(L) +
    bdiag(list(S1 - l1%*%t(l1),
               S2 - l2%*%t(l2),
               S3 - l3%*%t(l3),
               S4 - l4%*%t(l4),
               T5,
               T6))
  
  if(jac){
    return(implied_S[upper.tri(implied_S, diag = T)])
  }else{
    return(implied_S)
  }
}



get_loadings <- function(x, block_sizes) {
      #' Extract Loadings from Parameter Vector
      #'
      #' This function retrieves the loadings for each block of observed variables
      #' (manifest variables) based on the parameter vector `x` and the input block sizes.
      #'
      #' @param x A numeric vector of parameters. The first `y` values in `x` represent
      #'   the loadings to be extracted, where `y` is the total number of manifest
      #'   variables across all blocks.
      #' @param block_sizes A numeric vector where each element represents the number
      #'   of observed variables (manifest variables) in the corresponding block of
      #'   latent variables.
      #'
      #' @return A list of numeric vectors, where each vector contains the loadings
      #'   for the corresponding block of manifest variables.
      #'
      #' @examples
      #' block_sizes <- c(3, 2)
      #' x <- c(0.2, 0.4, 0.6, 0.1, 0.3)
      #' get_loadings(x, block_sizes)

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


get_correlation_coeff <- function(x, latent_variables, start_index) {
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
    #' latent_variables <- c("LV1", "LV2", "LV3")
    #' x <- c(0.5, 0.3, 0.2, 0.5, 0.4) # upper triangle values for a 3x3 matrix
    #' start_index <- 1
    #' get_correlation_coeff(x, latent_variables, start_index)

  # dim  matrix (number of exogeneous or endogeneous latent variables)
  dim <- length(latent_variables)
  #length upper values vector from this symetric matrix:
  length_correl <- dim * (dim - 1) / 2
  end_index <- start_index + length_correl - 1
  upper_values <- x[start_index:end_index]

  # Building of the correlation matrix
  P_matrix <- diag(1, dim, dim)
  if (dim > 1) {
    P_matrix[upper.tri(P_matrix)] <- upper_values
    P_matrix[lower.tri(P_matrix)] <- t(P_matrix)[lower.tri(P_matrix)]
  }

  return(P_matrix)
}


get_path_coeff <- function(x, list_linked_exo_endo, exo_or_endo_variable, initial_start_index) {
    #' Construct Path Coefficient Matrix from Parameter Vector
    #'
    #' This function generates a path coefficient matrix for endogenous and exogenous variables
    #' based on the parameter vector `x` and the provided relationships between variables.
    #'
    #' @param x A numeric vector of parameters. Contains path coefficient values to be extracted.
    #' @param list_linked_exo_endo A list where each element represents the linked exogenous or
    #'   endogenous variables for each endogenous variable. Keys correspond to variable names.
    #' @param exo_or_endo_variable A named numeric vector representing the exogenous or endogenous
    #'   variables. Used to define the structure of the path coefficient matrix.
    #' @param initial_start_index An integer specifying the index in `x` where the path coefficient
    #'   parameters start.
    #'
    #' @return A matrix where rows correspond to endogenous variables and columns to exogenous or
    #'   endogenous variables, with path coefficients filled in where applicable.
    #'
    #' @examples
    #' x <- c(0.8, 0.5, 0.3, 0.7, 0.8, 0.4)
    #' list_linked_exo_endo <- list(c(X1 = 1, X2 = 2), c(X3 = 3, X4 = 4))
    #' exo_or_endo_variable <- c(X1 = 1, X2 = 2, X3 = 3, X4 = 4)
    #' initial_start_index <- 1
    #' get_path_coeff(x, list_linked_exo_endo, exo_or_endo_variable, initial_start_index)

  # list of number of linked exos/endos for each endo variable
  # length_exo_endo <- sapply(list_linked_exo_endo, length)
  length_exo_endo <- sapply(list_linked_exo_endo,
                            function(sublist) {
                              ifelse((length(sublist) == 1 && sublist[[1]] == 0),
                                     return(0),
                                     return(length(sublist)))
                            })


  # Initialisation of index for the rows of the path coeff matrix
  start_index <- initial_start_index + c(0, cumsum(head(length_exo_endo, -1)))
  end_index <- start_index + length_exo_endo - 1

  # in the matrix, replace non zero value by x values
  Matrix_row <- mapply(
    function(linked_variable_i, start_ind, end_ind) {
      # Create a vector based on list_linked_exo_endo[i]
      result <- ifelse(names(exo_or_endo_variable) %in% names(linked_variable_i),
                       exo_or_endo_variable, 0)
      result[result != 0] <- x[start_ind:end_ind]
      return(result)
    },
    list_linked_exo_endo,                  # List of linked variables
    start_index,
    end_index,
    SIMPLIFY = FALSE
  )

  Matrix_path <- do.call(rbind, Matrix_row)

  return(Matrix_path)
}


build_formative_S_diag <- function(upper_values_composite_i, dim_block) {
    #' Construct Symmetric Formative S Matrix
    #'
    #' This function builds a symmetric covariance matrix (S_composite_i) for a formative block,
    #' using the upper triangular values provided in `upper_values_composite_i`.
    #'
    #' @param upper_values_composite_i A numeric vector containing the values for the upper
    #'   triangle (including the diagonal) of the symmetric matrix.
    #'
    #' @return A symmetric formative covariance matrix of dimensions `dim_block x dim_block`.
    #'
    #' @examples
    #' upper_values <- c(1, 0.3, 0.5, 1, 0.2, 1)
    #' dim_block <- 3
    #' build_formative_S_diag(upper_values, dim_block)

  # dim of the matrix extracted from the length of the upper values
  dim_block <- (-1 + sqrt(1 + 8 * length(upper_values_composite_i))) / 2
  S_composite_i <- matrix(0, dim_block, dim_block)
  S_composite_i[upper.tri(S_composite_i, diag = TRUE)] <- upper_values_composite_i
  S_composite_i[lower.tri(S_composite_i)] <- t(S_composite_i)[lower.tri(S_composite_i)]

  return(S_composite_i)
}


get_bdiag <- function(x, mode, block_sizes, loadings, initial_start_index_cov) {
      #' Construct Block Diagonal Covariance Structure
      #'
      #' This function generates a block diagonal matrix representing covariance structures
      #' for composite or formative blocks. Each diagonal block is computed based on the mode
      #' (formative or reflective) and the provided parameters.
      #'
      #' @param x A numeric vector of parameters from which covariance or variance values
      #'   are extracted.
      #' @param mode A vector of character strings indicating the mode for each block.
      #'   Possible values are "formative" or "reflective".
      #' @param block_sizes A numeric vector where each element represents the number
      #'   of emergent variables in the corresponding block.
      #' @param initial_start_index_cov An integer specifying the starting index in `x`
      #'   where covariance/variance parameters begin.
      #'
      #' @return A list of matrices, where each element represents a diagonal block of the
      #'   covariance structure. Blocks are either covariance matrices (formative blocks)
      #'   or diagonal matrices (reflective blocks).
      #'
      #' @examples
      #' x <- c(1, 0.3, 1, 0.2, 0.4, 1, 2, 0.5, 3)
      #' mode <- c("formative", "reflective")
      #' block_sizes <- c(2, 3)
      #' initial_start_index_cov <- 1
      #' get_bdiag(x, mode, block_sizes, initial_start_index_cov)

  J <- length(block_sizes)

  # list of lengths of the upper values in the cov matrix for the composite block i or
  # of the diagonal values for formative for each block
  lengths_values_cov <- ifelse(mode == "formative", (block_sizes^2 + block_sizes) / 2, block_sizes)

  start_indices_cov <- cumsum(c(initial_start_index_cov, head(lengths_values_cov, -1)))
  end_indices_cov <- start_indices_cov + lengths_values_cov - 1

  BDIAG <- list()

  for (i in seq_len(J)) {
    start_index_cov <- start_indices_cov[i]
    end_index_cov <- end_indices_cov[i]
    diag_or_upper_values_extracted <- x[start_index_cov:end_index_cov]

    if (mode[i] == "formative") {
      S_composite_i <- build_formative_S_diag(diag_or_upper_values_extracted)
      BDIAG [[i]] <- S_composite_i - loadings[[i]] %*% t(loadings[[i]])
    } else {
      BDIAG [[i]] <- diag(diag_or_upper_values_extracted)
    }
  }

  return(BDIAG)
}


get_bdiag_bis <- function(x, mode, block_sizes, initial_start_index_cov) {
    #' Construct Block Diagonal Covariance Structures for Formative and Reflective Blocks
    #'
    #' This function creates a block diagonal covariance matrix for both formative and reflective blocks
    #' using the parameter vector `x`, block sizes, and the starting index of covariance parameters.
    #'
    #' @param x A numeric vector of parameters from which covariance or variance values are extracted.
    #' @param mode A character vector indicating the mode ("formative" or "reflective") for each block.
    #' @param block_sizes A numeric vector where each element specifies the size (number of variables)
    #'   of each block.
    #' @param loadings A list of numeric vectors representing loadings for each block. Each vector corresponds
    #'   to the respective block size.
    #' @param initial_start_index_cov An integer indicating the starting index in `x` for covariance/variance
    #'   parameters.
    #'
    #' @return A list of block diagonal matrices, where each block corresponds to either:
    #'   - A covariance matrix for formative blocks.
    #'   - A diagonal matrix for reflective blocks.
    #'
    #' @examples
    #' x <- c(1, 0.3, 1, 0.2, 0.4, 1, 2, 0.5, 3)
    #' mode <- c("formative", "reflective")
    #' block_sizes <- c(2, 3)
    #' loadings <- list(c(0.8, 0.6), c(0.7, 0.5, 0.4))
    #' initial_start_index_cov <- 1
    #' get_bdiag_bis(x, mode, block_sizes, loadings, initial_start_index_cov)

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

  # Building list of formative matrices
  S_composites <- lapply(list_cov[mode == "formative"], build_formative_S_diag)

  # Building list of reflectives matrices
  reflective_blocks <- lapply(list_cov[mode == "reflective"], diag)

  BDIAG <- vector("list", J)
  BDIAG[mode == "formative"] <- S_composites
  BDIAG[mode == "reflective"] <- reflective_blocks

  return(BDIAG)

}



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




s_implied_bis <- function(x, block_sizes, mode, lengths_parameter, which_exo_endo, jac = TRUE){

  n <- which_exo_endo$ind_exo
  m <- which_exo_endo$ind_endo



  start_indices_in_x <- cumsum(c(1, head(lengths_parameter, -1)))

  ######################################################
  ####### mapping of the loadings from x ###############
  ######################################################

  loadings <- get_loadings(x, block_sizes = block_sizes)

  ##################################################################
  ####### mapping of the exogeneous correlation matrix from x ######
  ##################################################################

  P_EXO <- get_correlation_coeff(x,
                                 latent_variables = n,
                                 start_index = start_indices_in_x[2])

  ##################################################################
  ####### mapping of the path coefficients matrix G from x #########
  ##################################################################

  G <- get_path_coeff(x,
                      list_linked_exo_endo = which_exo_endo$Hi,
                      exo_or_endo_variable = n,
                      initial_start_index = start_indices_in_x[3])

  ##################################################################
  ####### mapping of the path coefficients matrix B from x #########
  ##################################################################

  B <- get_path_coeff(x,
                      list_linked_exo_endo = which_exo_endo$Ji,
                      exo_or_endo_variable = m,
                      initial_start_index = start_indices_in_x[4])

  ##################################################################
  ####### mapping of the endogeneous correlation matrix from x #####
  ##################################################################

  P_ENDO <- get_correlation_coeff(x,
                                  latent_variables = m,
                                  start_index = start_indices_in_x[5])

  ##################################################################
  ########################## Computation of PSI, R2 ################
  ##################################################################

  PSI <-  (diag(NROW(B)) - B)%*%P_ENDO%*%t((diag(NROW(B)) - B)) - G%*%P_EXO%*%t(G)
  R2 <- 1-diag(PSI)

  ########################################################################
  ########  Correlations between Latent/emergent variables   #############
  ########################################################################

  R <- rbind(cbind(P_EXO, P_EXO%*%t(G)%*%t(solve(diag(NROW(B)) - B))),
            cbind(solve(diag(NROW(B)) - B)%*%G%*%P_EXO, P_ENDO))

  ########################################################################
  ############## Compute the variance blocks #############################
  ########################################################################

  BDIAG <- get_bdiag_bis(x,
                         mode = mode,
                         block_sizes = block_sizes,
                         initial_start_index_cov = start_indices_in_x[6])

  # Get the residual variance for reflective blocks
  residual_variance <- lapply(BDIAG[mode == 'reflective'], diag)
  # Get the variance matrices for reflective blocks
  S_composites <- BDIAG[mode == 'formative']

  ########################################################################
  ############## Get omega ###############################################
  ########################################################################

  omega <- mapply(function(Sjj, lambda_j) solve(Sjj) %*% lambda_j,
                  S_composites, loadings[mode == "formative"], SIMPLIFY = FALSE)

  ########################################################################
  ####### Compute the implied covariance matrix implied by the model #####
  ########################################################################
  formative_loadings <- lapply(loadings[mode == "formative"], function(lambda) { lambda %*% t(lambda) })
  BDIAG[mode == 'formative'] <- mapply("-", S_composites, formative_loadings, SIMPLIFY = FALSE)

  L <- as.matrix(Matrix::bdiag(loadings))
  BDIAG <- as.matrix(Matrix::bdiag(BDIAG))

  implied_S <- L%*%R%*%t(L) + BDIAG


  out <- list(
    lambda = loadings,
    beta = B,
    gamma = G,
    psi = PSI,
    R2 = R2,
    residual_variance = residual_variance,
    S_composites = S_composites,
    omega = omega,
    P_EXO = P_EXO,
    P_ENDO = P_ENDO,
    R_LVM = R,
    SIGMA_IMPLIED  = implied_S
  )

  if(jac){
    return(implied_S[upper.tri(implied_S, diag = T)])
  }else{
    return(out)
  }

}
