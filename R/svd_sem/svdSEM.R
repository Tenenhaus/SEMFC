# #' Source Required R Files
# #'
# #' This function sources required R files for the svd-SEM implementation.
# #'
# #' @details
# #' Sources two R files:
# #' - lvm.R: Contains latent variable model functionality
# #' - d_LS.R: Contains least squares calculation utilities
# #'


source("R/svd_sem/correction.R")
source("R/svd_sem/lvm.R")
source("R/utils/d_LS.R")
source("R/utils/scale2.R")

#' @import Matrix
#' @importFrom readxl read_excel
#' Structural equation models with factors and composites with svd-SEM
#' @param A  A list that contains the \eqn{J} blocks of variables \eqn{X_1, X_2, ..., X_J}.
#' @param scale  If scale = TRUE, each block is standardized to zero means and unit variances.
#' @param bias A logical value for biaised or unbiaised estimator of the var/cov (default: bias = TRUE).
#' @param mode a vector of lenght J indication the mode for each block formative or reflective.
#' @return \item{a}{A list of \eqn{J} elements. Each element of \eqn{a} is a matrix that contains the first singular vector for each block.}
#' @references Tenenhaus M., Tenenhaus A. and Groenen PJF (2017), Regularized generalized canonical correlation analysis: A framework for sequential multiblock component methods, Psychometrika, in press
#' @title Structural Equation Modeling with Factors and Composites (svdSEM) 
#' @examples
#' #############
#' # Example 1 #
#' #############
#' ECSI = as.data.frame(read_excel("mobil.xls"))/10
#' A = list(CUSTOMER_E = ECSI[, c("CUEX1", "CUEX2", "CUEX3")],
#'          PERC_QUAL  = ECSI[, c("PERQ1", "PERQ2", "PERQ3", "PERQ4", 
#'                                "PERQ5", "PERQ6", "PERQ7")],
#'          PERC_VALUE = ECSI[, c("PERV1", "PERV2")],
#'          CUSTOMER_S = ECSI[, c("CUSA1", "CUSA2", "CUSA3")],
#'          CUSTOMER_L = ECSI[, c("CUSL1", "CUSL2", "CUSL3")])
#'          
#' C = matrix(c(0, 0, 0, 0, 0,
#'              1, 0, 0, 0, 0,
#'              1, 1, 0, 0, 0,
#'              1, 1, 1, 0, 0,
#'              0, 0, 0, 1, 0), 5, 5, byrow = FALSE)
#' 
#' fit = svdSEM(A, C, scale = FALSE, 
#'              mode = rep("reflective", length(A)),
#'             bias = FALSE)
#'             
#' @export rgcca
svdSEM <- function(A, C, scale = TRUE,
                   mode = rep("formative", length(A)), 
                   bias = FALSE, pen=FALSE){

    # vecteur qui code le nombre des variables observes de chaque bloc
    pjs <- sapply(A, NCOL)
    # nombre d'individu
    nb_row <- NROW(A[[1]])
    
    #-------------------------------------------------------
    blocks = A
    if(scale){
      A = lapply(A, function(x) scale2(x, bias = bias))
    }else{
      A = lapply(A, function(x) scale2(x, bias = bias, scale = FALSE))
    }
    

    lambda = list()
    omega = list()
    nb_ind <- NROW(A[[1]])
    J <- length(A)



    library(PMA)


    psvd <- function(x){
      data <- t(A[[x]])%*%Reduce("cbind", A[-x])%*%t(t(A[[x]])%*%Reduce("cbind", A[-x]))
      # data <- t(A[[x]])%*%Reduce("cbind", A[-x])  u et v peuvent ne pas etre les memes attention



      cv.out <- SPC.cv(data, sumabsvs=seq(1, sqrt(ncol(data)), len=20))

      out <- SPC(data, sumabsv=cv.out$bestsumabsv, K=1, v=cv.out$v.init)

      return(out$v)


    }

    psvd2 <- function(x){

      data <- t(A[[x]])%*%Reduce("cbind", A[-x])



      cv.out <- PMD.cv(data, type = 'standard',sumabss=seq(1/sqrt(nrow(data)), 1, len=50))

      out <- PMD(data, type = 'standard', sumabs=cv.out$bestsumabs, K=1, v=cv.out$v.init)

      return(out$u)


    }

    if (pen==FALSE){

      # Extract first singular vector for each block.
      a = sapply(1:J,
               function(x)
                 svd(t(A[[x]])%*%Reduce("cbind", A[-x]),
                     nu = 1, nv = 1)$u,
               simplify = FALSE
               )
    } else {
      a = sapply(1:J,
               function(x)
                 psvd(x),
               simplify = FALSE
               )

      a2 = sapply(1:J,
                   function(x)
                     psvd2(x),
                   simplify = FALSE
                   )

      asvd = sapply(1:J,
                   function(x)
                     svd(t(A[[x]])%*%Reduce("cbind", A[-x]),
                         nu = 1, nv = 1)$u,
                   simplify = FALSE
                   )

    }





    
    #check for sign inversion 
    a = lapply(a, function(x) {if (x[1]>0) {x<-x} else {x<--x}})
    names(a) = names(A)

    #Compute disattenuation 
    d <- rep(0, J)
    
    for(j in seq_len(J)){
      d[j] <- correction(A[[j]], a[[j]], 
                         mode = mode[j], 
                         bias = bias)$d
    }
    
    names(d) = names(A)
    
    #... and apply the correction
    lambda = mapply("*", a, d, SIMPLIFY = FALSE)
    
    for (j in which(mode == "formative")){
        omega[[j]] = solve(cov2(A[[j]], bias = bias))%*%a[[j]]*d[j]
    }
    
    if(any(mode == "formative")) 
      names(omega) = names(A)[which(mode=="formative")]
    
        
    for (b in seq_len(J))
      rownames(a[[b]]) = rownames(lambda[[b]]) = colnames(A[[b]])  
      

    # reliability coefficients 
    reliability_coef = d^2/mapply("%*%", 
                                  lapply(a, t), 
                                  mapply("%*%", 
                                  lapply(A, function(x) cov2(x, bias= bias)), 
                                         a, SIMPLIFY = FALSE)
                                  )
        
      
    # rho_jk (Expected value : -1 <=rho_jk <=1)
    Phat = matrix(0, J, J)
    for (j in seq_len(J)){
      for (k in j:J){
            Phat[j, k] = t(a[[j]])%*%cov2(A[[j]], A[[k]], bias = bias)%*%a[[k]]/(d[j]*d[k])
      }
    }
    
    Phat = Phat + t(Phat)
    diag(Phat) = 1
    colnames(Phat) = rownames(Phat) = names(A)
    
    #Call lvm() for the structural model
    lv = lvm(Phat, C)
    
    var_MVs = lapply(A, function(x) diag(cov2(x, bias = bias)))
    
    #standardized loadings     
    # Expected value : -1 <= cor(y_jh, eta_j) <=1
    std_lambda = mapply("/", lambda, lapply(var_MVs, sqrt),  SIMPLIFY = FALSE)
    
    #residual variance    
    residual_variance = mapply("-", var_MVs, lapply(lambda, function(x) x^2), 
                               SIMPLIFY = FALSE)
    
    #Compute the implied covariance matrix implied by the model
    BDIAG = list()

    for (i in seq_len(J)){
      if(mode[i] == "formative"){
        BDIAG [[i]] = cov2(A[[i]], bias = bias)-lambda[[i]]%*%t(lambda[[i]])
      }else{
        BDIAG [[i]] = diag(drop(residual_variance[[i]]))
      }
    }
    
    LAMBDA = as.matrix(Matrix::bdiag(lambda))
    BDIAG = as.matrix(Matrix::bdiag(BDIAG))
    SIGMA_LVM = LAMBDA%*%lv$R_LVM%*%t(LAMBDA) + BDIAG
    
    var_Names = colnames(Reduce("cbind", A))
    dimnames(SIGMA_LVM) = list(var_Names, var_Names) 
    
    #Measure of goodness-of(fit)
    T_LS = d_LS(cov2(Reduce("cbind", A), bias = bias), SIGMA_LVM)
    

    out <- list(a = a, 
                lambda = lambda,
                std_lambda = std_lambda,
                omega = omega,
                gr = lv$gr,
                beta = lv$BETA,
                gamma = lv$GAMMA,
                R2 = lv$R2,
                psi = lv$PSI,
                d = d,
                Ptilde = Phat,
                P_EXO = lv$P_EXO,
                P_ENDO = lv$P_ENDO,
                P_IMPLIED = lv$R_LVM,
                SIGMA_IMPLIED = SIGMA_LVM,
                T_LS = T_LS, 
                reliability_coef = reliability_coef,
                residual_variance = residual_variance[mode == 'reflective'],
                blocks = blocks,
                mode = mode,
                bias = bias,
                scale = scale,
                C = C)
        
    return(out)
}
    
    