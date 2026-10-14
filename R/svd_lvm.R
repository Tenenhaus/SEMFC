
#' Latent Variable Model (LVM) Estimation
#'
#' Estimates the parameters of a latent variable model (LVM) using the provided
#' correlation matrix  and the adjacency matrix of the structural model.
#'
#' @param R Numeric matrix representing the correlation matrix.
#' @param C Numeric matrix representing the adjacency matrix of the structural model.
#'
#' @return A list containing:
#'   \item{gr}{Graph object representing the structural model.}
#'   \item{BETA}{Matrix of path coefficients between endogenous variables.}
#'   \item{GAMMA}{Matrix of path coefficients from exogenous to endogenous variables.}
#'   \item{PSI}{Residual covariance matrix of endogenous variables.}
#'   \item{R2}{Vector of R-squared values for endogenous variables.}
#'   \item{P_EXO}{Implied covariance matrix of exogenous latent variables.}
#'   \item{P_ENDO}{Implied covariance matrix of endogenous latent variables.}
#'   \item{R_LVM}{Implied covariance matrix of all latent variables.}
#'
#' @details
#' The function performs the following steps:
#' - Identifies exogenous and endogenous latent variables using the adjacency matrix.
#' - Computes path coefficients (BETA, GAMMA) and residual covariance matrix (PSI).
#' - Constructs the implied covariance matrix of latent variables (R_LVM).
#' - Handles both directed acyclic graphs (DAGs) and non-DAGs.
#'
#' @importFrom igraph graph_from_adjacency_matrix is_dag
#'
#' @keywords internal
lvm <- function(R, C){
  gr = igraph::graph_from_adjacency_matrix(C)
  which_exo_endo = ind_exo_endo(C)
  rownames(C) = colnames(C) = rownames(R)
  xx = list()
  xy = list()
  bg = list()
  b = list()
  g = list()
  tC = t(C)
  H = which_exo_endo$ind_exo
  J = which_exo_endo$ind_endo
  
  BETA = as.matrix(tC[J, J, drop = F])
  GAMMA = as.matrix(tC[J, H, drop = F])
  PSI = matrix(0, NCOL(BETA), NCOL(BETA))

  for (i in 1:length(which_exo_endo$Ji)){
    # les endo lie a l'endo i
    Ji = which_exo_endo$Ji[[i]]
    # les exos lies a l'endo i
    Hi = which_exo_endo$Hi[[i]]
    l_j = length(Ji)
    l_h = length(Hi)
    l = length(c(Hi, Ji))
    
    config1 = any(Ji != 0) & any(Hi != 0)
    config2 = any(Ji != 0) & any(Hi == 0)
    config3 = any(Ji == 0) & any(Hi != 0)
    config4 = any(Ji == 0) & any(Hi == 0)
    
    if(!igraph::is_dag(gr) & config1){
      xx = matrix(NA, l, l)    
      xx[1:l_j, 1:l_j] = R[Ji, H]%*%solve(R[H, H])%*%R[H, Ji]
      xx[1:l_j, (l_j+1):l] = R[Ji, Hi]
      xx[(l_j+1):l, 1:l_j] = R[Hi, Ji]
      xx[(l_j+1):l, (l_j+1):l] = R[Hi, Hi]
      
      xy = c(R[Ji, H]%*%solve(R[H, H])%*%R[H, J[[i]]], R[Hi, J[[i]]])
      bg[[i]] = solve(xx)%*%xy
    }
    
    if(!igraph::is_dag(gr) & config2){
      xx = R[Ji, H]%*%solve(R[H, H])%*%R[H, Ji]
      xy = R[Ji, H]%*%solve(R[H, H])%*%R[H, J[[i]]]
      bg[[i]] = solve(xx)%*%xy
      
    }
    
    if(!igraph::is_dag(gr) & config3){
      xx = R[Hi, Hi]
      xy = R[Hi, J[[i]]]
      bg[[i]] = solve(xx)%*%xy
    }
    
    if(!igraph::is_dag(gr) & config4){
      print(paste("the ", i, "th structural model is not properly specified", 
                  sep = "")
            )
      break 
    }
    
    if(igraph::is_dag(gr) & config1){
    #   cov vect(endo-i,exo-i)
    xx = matrix(NA, l, l)
    xx[1:l_j, 1:l_j] = R[Ji, Ji]
    xx[1:l_j, (l_j+1):l] = R[Ji, Hi]
    xx[(l_j+1):l, 1:l_j] = R[Hi, Ji]
    xx[(l_j+1):l, (l_j+1):l] = R[Hi, Hi] 
    xy = c(R[Ji, J[i]], R[Hi, J[i]])
    bg[[i]] = solve(xx)%*%xy
    }
    
    if(igraph::is_dag(gr) & config2){
      
      xx = R[Ji, Ji]
      xy = R[Ji, J[i]]
      bg[[i]] = solve(xx)%*%xy
    }
    
    if(igraph::is_dag(gr) & config3){
      xx = R[Hi, Hi] 
      xy = R[Hi, J[i]]
      bg[[i]] = solve(xx)%*%xy
    }
    

    if(igraph::is_dag(gr) & config4){
      print(paste("the ", i, "th structural model is not properly specified", 
                  sep = "")
      )
      break 
    }

    ifelse(Ji == 0, 
           yes = {b[[i]] = 0},
           no = {b[[i]] = list(bg[[i]][1:length(Ji)])} 
    )
    
    ifelse(Hi == 0, 
      yes = {g[[i]] = 0},
      no = {ifelse(Ji == 0,
            yes = {g[[i]] = list(bg[[i]])},
            no = {g[[i]] = list(bg[[i]][seq(length(bg[[i]]))[-seq(length(Ji))]])})
           })
    
    BETA[i, ] = replace(BETA[i, ], BETA[i, ] == 1, b[[i]][[1]])
    GAMMA[i, ] = replace(GAMMA[i, ], GAMMA[i , ] == 1, g[[i]][[1]])
    

  }

  PI = solve(diag(NROW(BETA))- BETA)
  if(!igraph::is_dag(gr)){
    PSI = (diag(NROW(BETA))- BETA)%*%R[J,J]%*%t(diag(NROW(BETA))- BETA) -
      GAMMA%*%R[H,H]%*%t(GAMMA)
  } else {
    D <- diag(diag(ncol(BETA)) - (PI%*%GAMMA%*%R[H, H]%*%t(GAMMA)%*%t(PI)))
    diag_PSI <- drop(solve(PI*PI)%*%D)
    PSI <- if(length(diag_PSI) == 1) {
      matrix(diag_PSI, 1, 1)
    } else {
      diag(diag_PSI)
    }

  }
  dimnames(PSI) <- list(
    paste0(".", rownames(BETA)),
    paste0(".", rownames(BETA))
  )

  R2 = 1-diag(PSI)

  R_LVM = matrix(0, NCOL(C), NCOL(C))
  
  if(!igraph::is_dag(gr)){
    R_LVM[H, H] = R[H, H]
    R_LVM[J, J] = R[J, J]
    R_LVM[H, J] = R[H, H]%*%t(GAMMA)%*%t(PI)
    R_LVM[J, H] = PI%*%GAMMA%*%R[H, H]
  }else{
    R_LVM[H, H] = R[H, H]
    R_LVM[H, J] = R[H, H]%*%t(GAMMA)%*%t(PI)
    R_LVM[J, H] = PI%*%GAMMA%*%R[H, H]
    R_LVM[J, J] = PI%*%(GAMMA%*%R[H, H]%*%t(GAMMA) + PSI)%*%t(PI)
  }

  dimnames(R_LVM) <- list(rownames(C), rownames(C))
  

  return(list(gr = gr, 
              BETA = BETA, GAMMA = GAMMA, 
              PSI = PSI, 
              R2 = R2,
              P_EXO = R_LVM[H, H],
              P_ENDO = R_LVM[J, J],
              R_LVM = R_LVM))
}




#' Structure information for exogenous and endogenous latent variables
#'
#' @param C matrix filled with 1 and 0 with the same
#' number of columns and rows and without any 1 on the diagonal
#'
#' @return list of H, J, index of the endogenous variables and
#' index of the exogenous variables of the C matrix
#'
#'
#' @examples
#' \dontrun{
#' C = matrix(c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,1,
#' 0,0,1,1,1,0),6,6)
#' ind_exo_endo(C)
#' }
#' @keywords internal



ind_exo_endo = function(C){
  #Creation of copy
  DBIS = DBIS2 = C
  ind_exo = which(colSums(C) == 0)
  #We look at which columns are only composed of 0 and we transform
  #the column corresponding to 0
  DBIS[ind_exo, ] = 0
  #We create a new matrix without any exogenous variables
  omega = DBIS[ , -ind_exo]
  omega = as.matrix(omega)
  #We display the J
  Ji = sapply(1:length(which(colSums(C) != 0)),
              function(x) which(omega[ , x] != 0),
              simplify = FALSE
  )

  #We transform all the rows whose columns are different from 0 into 0
  #to end up with a matrix where only the exogenous variables have ones
  ind_endo = which(colSums(C) != 0)

  DBIS2[ind_endo,] = 0
  gamma = DBIS2[ , -ind_exo]
  gamma = as.matrix(gamma)

  #We display H
  Hi = sapply(1:length(which(colSums(C) != 0)),
              function(x) which(gamma[ , x] != 0),
              simplify = FALSE
  )

  Hi[sapply(Hi, function(x) length(x) == 0)] <- 0
  Ji[sapply(Ji, function(x) length(x) == 0)] <- 0

  return(list(Hi = Hi, Ji = Ji,ind_endo = ind_endo,ind_exo = ind_exo))
}


