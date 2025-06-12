#' ind_exo_endo
#'
#' @param C matrix filled with 1 and 0 with the same 
#' number of columns and rows and without any 1 on the diagonal
#'
#' @return list of H, J, index of the endogenous variables and 
#' index of the exogenous variables of the C matrix
#' 
#' @export
#'
#' @examples
#' C = matrix(c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,1,
#' 0,0,1,1,1,0),6,6)
#' ind_exo_endo(C)




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
