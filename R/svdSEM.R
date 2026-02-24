

#' @import Matrix
#' @param A  A list that contains the \eqn{J} blocks of variables \eqn{X_1, X_2, ..., X_J}.
#' @param C  An adjacency matrix representing the structural model among latent variables.
#' @param scale  If scale = TRUE, each block is standardized to zero means and unit variances.
#' @param bias A logical value for biaised or unbiaised estimator of the var/cov (default: bias = TRUE).
#' @param mode a vector of lenght J indication the mode for each block formative or reflective.
#' @return A list containing the following elements:
#' \describe{
#'   \item{a}{A list of J elements. Each element contains the first singular vector for each block.}
#'   \item{lambda}{A list of loadings for each block.}
#'   \item{std_lambda}{A list of standardized loadings for each block.}
#'   \item{std_omega}{A list of standardized weights for formative blocks only.}
#'   \item{omega}{A list of weights for formative blocks only.}
#'   \item{gr}{A directed graph representing the structural model.}
#'   \item{beta}{A matrix of structural coefficients among endogenous latent variables.}
#'   \item{gamma}{A matrix of structural coefficients from exogenous to endogenous latent variables.}
#'   \item{R2}{A vector of R-squared coefficients for endogenous latent variables.}
#'   \item{psi}{A variance-covariance matrix of structural errors.}
#'   \item{d}{A vector of disattenuation factors for each block.}
#'   \item{Ptilde}{An estimated correlation matrix among latent variables.}
#'   \item{P_EXO}{A correlation matrix of exogenous latent variables.}
#'   \item{P_ENDO}{A correlation matrix of endogenous latent variables.}
#'   \item{P_IMPLIED}{A model-implied correlation matrix for latent variables.}
#'   \item{SIGMA_IMPLIED}{A model-implied covariance matrix for observed variables.}
#'   \item{T_LS}{A model fit measure (least squares discrepancy).}
#'   \item{reliability_coef}{A vector of reliability coefficients for each block.}
#'   \item{residual_variance}{A list of residual variances for reflective blocks.}
#'   \item{blocks}{The input data (list of blocks).}
#'   \item{mode}{A vector of measurement modes for each block.}
#'   \item{bias}{A logical value indicating whether biased estimation was used.}
#'   \item{scale}{A logical value indicating whether data were standardized.}
#'   \item{C}{An adjacency matrix of the structural model.}
#' }
#' @references Tenenhaus M., Tenenhaus A. and Groenen PJF (2017), Regularized generalized canonical correlation analysis: A framework for sequential multiblock component methods, Psychometrika, in press
#' @title Structural Equation Modeling with Factors and Composites (svdSEM) 
#' @examples
#' #############
#' # Example 1 #
#' #############
#' data(ECSI)
#' ECSI = ECSI/10
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
#' @export
svdSEM <- function(A, C, scale = TRUE,
                   mode = rep("formative", length(A)), 
                   bias = FALSE){

    #-------------------------------------------------------
    blocks <- A
    if(scale){
      A <- lapply(A, function(x) scale2(x, bias = bias))
    }else{
      A <- lapply(A, function(x) scale2(x, bias = bias, scale = FALSE))
    }
    

    lambda <- list()
    omega <- list()
    std_omega <- list()
    nb_ind <- NROW(A[[1]])
    J <- length(A)


    # Extract first singular vector for each block.
    a <- sapply(1:J,
               function(x) 
                 svd(t(A[[x]])%*%Reduce("cbind", A[-x]), 
                     nu = 1, nv = 1)$u, 
               simplify = FALSE
               )
    
    #check for sign inversion 
    a <- lapply(a, function(x) {if (x[1]>0) {x<-x} else {x<--x}})
    names(a) <- names(A)

    #Compute disattenuation 
    d <- rep(0, J)
    
    for(j in seq_len(J)){
      d[j] <- correction(A[[j]], a[[j]], 
                         mode = mode[j], 
                         bias = bias)$d
    }
    
    names(d) <- names(A)
    
    #... and apply the correction
    lambda <- mapply("*", a, d, SIMPLIFY = FALSE)
    
    for (j in which(mode == "formative")){
        omega[[j]] <- solve(cov2(A[[j]], bias = bias))%*%a[[j]]*d[j]
    }
    omega <- omega[!sapply(omega, is.null)]

    
    if(any(mode == "formative")) 
      names(omega) <- names(A)[which(mode=="formative")]
    
        
    for (b in seq_len(J))
      rownames(a[[b]]) <- rownames(lambda[[b]]) <- colnames(A[[b]])
      

    # reliability coefficients 
    reliability_coef <- d^2/mapply("%*%",
                                  lapply(a, t), 
                                  mapply("%*%", 
                                  lapply(A, function(x) cov2(x, bias= bias)), 
                                         a, SIMPLIFY = FALSE)
                                  )
        
      
    # rho_jk (Expected value : -1 <=rho_jk <=1)
    Phat <- matrix(0, J, J)
    for (j in seq_len(J)){
      for (k in j:J){
            Phat[j, k] <- t(a[[j]])%*%cov2(A[[j]], A[[k]], bias = bias)%*%a[[k]]/(d[j]*d[k])
      }
    }
    
    Phat <- Phat + t(Phat)
    diag(Phat) <- 1
    colnames(Phat) <- rownames(Phat) <- names(A)
    
    #Call lvm() for the structural model
    lv <- lvm(Phat, C)
    
    var_MVs <- lapply(A, function(x) diag(cov2(x, bias = bias)))
    
    #standardized loadings     
    # Expected value : -1 <= cor(y_jh, eta_j) <=1
    std_lambda <- mapply("/", lambda, lapply(var_MVs, sqrt),  SIMPLIFY = FALSE)


    for (j in which(mode == "formative")){
        std_omega[[j]] <- solve(cov2(A[[j]], bias = bias))%*%std_lambda[[j]]
    }
    std_omega <- std_omega[!sapply(std_omega, is.null)]

    
    #residual variance    
    residual_variance <- mapply("-", var_MVs, lapply(lambda, function(x) x^2),
                               SIMPLIFY = FALSE)
    
    #Compute the implied covariance matrix implied by the model
    BDIAG <- list()

    for (i in seq_len(J)){
      if(mode[i] == "formative"){
        BDIAG [[i]] <- cov2(A[[i]], bias = bias)-lambda[[i]]%*%t(lambda[[i]])
      }else{
        BDIAG [[i]] <- diag(drop(residual_variance[[i]]))
      }
    }
    
    LAMBDA <- as.matrix(Matrix::bdiag(lambda))
    BDIAG <- as.matrix(Matrix::bdiag(BDIAG))
    SIGMA_LVM <- LAMBDA%*%lv$R_LVM%*%t(LAMBDA) + BDIAG
    
    var_Names <- colnames(Reduce("cbind", A))
    dimnames(SIGMA_LVM) <- list(var_Names, var_Names)
    
    #Measure of goodness-of(fit)
    T_LS <- d_LS(cov2(Reduce("cbind", A), bias = bias), SIGMA_LVM)

    lambda <- lapply(lambda, function(x) setNames(as.vector(x), rownames(x)))
    residual_variance <- lapply(residual_variance, function(x) setNames(as.vector(x), paste0(".", rownames(x))))


    out <- list(a = a, 
                lambda = lambda,
                std_lambda = std_lambda,
                omega = omega,
                std_omega = std_omega,
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
    
    