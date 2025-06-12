# print.rgcca <- function(x,verbose=FALSE,...) {
#   nbloc <- length(x$Y)
#   cat("Outputs for RGCCA: \n")
# }

#' Regularized Generalized Canonical Correlation Analysis (RGCCA) is a generalization
#' of regularized canonical correlation analysis to three or more sets of variables. 
#' Given \eqn{J} matrices \eqn{X_1, X_2, ..., X_J} that represent 
#' \eqn{J} sets of variables observed on the same set of \eqn{n} individuals. The matrices 
#' \eqn{X_1, X_2, ..., X_J} must have the same number of rows, 
#' but may (and usually will) have different numbers of columns. The aim of RGCCA is to study 
#' the relationships between these \eqn{J} blocks of variables. It constitutes a general 
#' framework for many multi-block data analysis methods. It combines the power of 
#' multi-block data analysis methods (maximization of well identified criteria) 
#' and the flexibility of PLS path modeling (the researcher decides which blocks 
#' are connected and which are not). Hence, the use of RGCCA requires the construction 
#' (user specified) of a design matrix \eqn{C}, that characterize 
#' the connections between blocks. Elements of the symmetric design matrix \eqn{C = (c_{jk})} 
#' is equal to 1 if block \eqn{j} and block \eqn{k} are connected, and 0 otherwise.
#' The function rgcca() implements a monotonically convergent algorithm (i.e. the bounded 
#' criteria to be maximized increases at each step of the iterative procedure) that is very 
#' similar to the PLS algorithm proposed by Herman Wold and finds at convergence a stationnary point 
#' of the RGCCA optimization problem. . Moreover, depending on the 
#' dimensionality of each block \eqn{X_j}, \eqn{j = 1, ..., J}, the primal (when \eqn{n > p_j}) algorithm or 
#' the dual (when \eqn{n < p_j}) algorithm is used (see Tenenhaus et al. 2015). 
#' Moreover, by deflation strategy, rgcca() allow to compute several RGCCA block 
#' components (specified by ncomp) for each block. Within each block, block components are guaranteed to 
#' be orthogonal using the deflation procedure. The so-called symmetric deflation is considered in
#' this implementation, i.e. each block is deflated with respect to its own component(s).
#' It should be noted that the numbers of components per block can differ from one block to another. 
#' @param A  A list that contains the \eqn{J} blocks of variables \eqn{X_1, X_2, ..., X_J}.
#' @param scale  If scale = TRUE, each block is standardized to zero means and unit variances and then divided by the square root of its number of variables (default: TRUE).
#' @param verbose  If verbose = TRUE, the progress will be report while computing (default: TRUE).
#' @param init The mode of initialization to use in RGCCA algorithm. The alternatives are either by Singular Value Decompostion ("svd") or random ("random") (Default: "svd").
#' @param bias A logical value for biaised or unbiaised estimator of the var/cov (default: bias = TRUE).
#' @param tol The stopping value for convergence.
#' @return \item{Y}{A list of \eqn{J} elements. Each element of \eqn{Y} is a matrix that contains the RGCCA components for the corresponding block.}
#' @return \item{a}{A list of \eqn{J} elements. Each element of \eqn{a} is a matrix that contains the outer weight vectors for each block.}
#' @return \item{astar}{A list of \eqn{J} elements. Each element of astar is a matrix defined as Y[[j]][, h] = A[[j]]\%*\%astar[[j]][, h].}
#' @return \item{crit}{A vector that contains the values of the criteria across iterations.}
#' @return \item{primal_dual}{A \eqn{1 * J} vector that contains the formulation ("primal" or "dual") applied to each of the \eqn{J} blocks within the RGCCA alogrithm} 
#' @return \item{AVE}{indicators of model quality based on the Average Variance Explained (AVE): AVE(for one block), AVE(outer model), AVE(inner model).}
#' @references Tenenhaus M., Tenenhaus A. and Groenen PJF (2017), Regularized generalized canonical correlation analysis: A framework for sequential multiblock component methods, Psychometrika, in press
#' @references Tenenhaus A., Philippe C., & Frouin V. (2015). Kernel Generalized Canonical Correlation Analysis. Computational Statistics and Data Analysis, 90, 114-131.
#' @references Tenenhaus A. and Tenenhaus M., (2011), Regularized Generalized Canonical Correlation Analysis, Psychometrika, Vol. 76, Nr 2, pp 257-284.
#' @references Schafer J. and Strimmer K., (2005), A shrinkage approach to large-scale covariance matrix estimation and implications for functional genomics. Statist. Appl. Genet. Mol. Biol. 4:32.
#' @title Regularized Generalized Canonical Correlation Analysis (RGCCA) 
#' @examples
#' #############
#' # Example 1 #
#' #############
#' data(Russett)
#' X_agric =as.matrix(Russett[,c("gini","farm","rent")])
#' X_ind = as.matrix(Russett[,c("gnpr","labo")])
#' X_polit = as.matrix(Russett[ , c("demostab", "dictator")])
#' A = list(X_agric, X_ind, X_polit)
#' #Define the design matrix (output = C) 
#' C = matrix(c(0, 0, 1, 0, 0, 1, 1, 1, 0), 3, 3)
#' result.rgcca = rgcca(A, C, tau = c(1, 1, 1), scheme = "factorial", scale = TRUE)
#' lab = as.vector(apply(Russett[, 9:11], 1, which.max))
#' plot(result.rgcca$Y[[1]], result.rgcca$Y[[2]], col = "white", 
#'      xlab = "Y1 (Agric. inequality)", ylab = "Y2 (Industrial Development)")
#' text(result.rgcca$Y[[1]], result.rgcca$Y[[2]], rownames(Russett), col = lab, cex = .7)
#'
#' #############
#' # Example 2 #
#' #############
#' data(Russett)
#' X_agric =as.matrix(Russett[,c("gini","farm","rent")])
#' X_ind = as.matrix(Russett[,c("gnpr","labo")])
#' X_polit = as.matrix(Russett[ , c("inst", "ecks", "death", 
#'                                  "demostab", "dictator")])
#' A = list(X_agric, X_ind, X_polit, cbind(X_agric, X_ind, X_polit))
#' 
#' #Define the design matrix (output = C) 
#' C = matrix(c(0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0), 4, 4)
#' result.rgcca = rgcca(A, C, tau = c(1, 1, 1, 0), ncomp = rep(2, 4), 
#'                      scheme = function(x) x^4, scale = TRUE) # HPCA
#' lab = as.vector(apply(Russett[, 9:11], 1, which.max))
#' plot(result.rgcca$Y[[4]][, 1], result.rgcca$Y[[4]][, 2], col = "white", 
#'      xlab = "Global Component 1", ylab = "Global Component 2")
#' text(result.rgcca$Y[[4]][, 1], result.rgcca$Y[[4]][, 2], rownames(Russett), 
#'      col = lab, cex = .7)
#'  
#' \dontrun{
#' ######################################
#' # example 3: RGCCA and leave one out #
#' ######################################
#' Ytest = matrix(0, 47, 3)
#' X_agric =as.matrix(Russett[,c("gini","farm","rent")])
#' X_ind = as.matrix(Russett[,c("gnpr","labo")])
#' X_polit = as.matrix(Russett[ , c("demostab", "dictator")])
#' A = list(X_agric, X_ind, X_polit)
#' #Define the design matrix (output = C) 
#' C = matrix(c(0, 0, 1, 0, 0, 1, 1, 1, 0), 3, 3)
#' result.rgcca = rgcca(A, C, tau = rep(1, 3), ncomp = rep(1, 3), 
#'                      scheme = "factorial", verbose = TRUE)
#'                      
#' for (i in 1:nrow(Russett)){
#'  B = lapply(A, function(x) x[-i, ])
#'  B = lapply(B, scale2)
#'  resB = rgcca(B, C, tau = rep(1, 3), scheme = "factorial", scale = FALSE, verbose = FALSE)
#'  #  look for potential conflicting sign among components within the loo loop.
#'  for (k in 1:length(B)){
#'    if (cor(result.rgcca$a[[k]], resB$a[[k]]) >= 0) 
#'      resB$a[[k]] = resB$a[[k]] else resB$a[[k]] = -resB$a[[k]]
#'  }
#'  Btest =lapply(A, function(x) x[i, ])
#'  Btest[[1]]=(Btest[[1]]-attr(B[[1]],"scaled:center")) /
#'                  (attr(B[[1]],"scaled:scale"))/sqrt(NCOL(B[[1]]))
#'  Btest[[2]]=(Btest[[2]]-attr(B[[2]],"scaled:center")) / 
#'                  (attr(B[[2]],"scaled:scale"))/sqrt(NCOL(B[[2]]))
#'  Btest[[3]]=(Btest[[3]]-attr(B[[3]],"scaled:center")) / 
#'                  (attr(B[[3]],"scaled:scale"))/sqrt(NCOL(B[[3]]))
#'  Ytest[i, 1] = Btest[[1]]%*%resB$a[[1]]
#'  Ytest[i, 2] = Btest[[2]]%*%resB$a[[2]]
#'  Ytest[i, 3] = Btest[[3]]%*%resB$a[[3]]
#' }
#' lab = apply(Russett[, 9:11], 1, which.max)
#' plot(result.rgcca$Y[[1]], result.rgcca$Y[[2]], col = "white", 
#'      xlab = "Y1 (Agric. inequality)", ylab = "Y2 (Ind. Development)")
#' text(result.rgcca$Y[[1]], result.rgcca$Y[[2]], rownames(Russett), 
#'      col = lab, cex = .7)
#' text(Ytest[, 1], Ytest[, 2], substr(rownames(Russett), 1, 1), 
#'      col = lab, cex = .7)
#'}
#' @export rgcca


rgccac <- function(A, C, scale = TRUE , 
                   mode = rep("formative", length(A)), 
                   init="svd", bias = FALSE, tol = 1e-8, verbose=TRUE){
    pjs <- sapply(A, NCOL)
    nb_row <- NROW(A[[1]])
    
    #-------------------------------------------------------
    
    if(scale){
      A = lapply(A, function(x) scale2(x, bias = bias))
    }else{
      A = lapply(A, function(x) scale2(x, bias = bias, scale = FALSE))
    }
    
    AVE_X = list()
    AVE_outer <- vector()
    lambda = list()
    nb_ind <- NROW(A[[1]])
    J <- length(A)
    improper_solution = list()
    
    result <- rgccak(A, C = 1-diag(J), tau = rep(1, J), scheme = "factorial", 
                     init = init, bias = bias, tol = tol, verbose=verbose)
    
    Y <- NULL
    for (b in 1:J) Y[[b]] <- result$Y[, b, drop = FALSE]
      
    #Average Variance Explained (AVE) per block
    for (j in 1:J) AVE_X[[j]] =  mean(cor(A[[j]], Y[[j]])^2)
        
    #AVE outer 
    AVE_outer <- sum(pjs * unlist(AVE_X))/sum(pjs)
                
    AVE <- list(AVE_X = AVE_X, 
                AVE_outer = AVE_outer,
                AVE_inner = result$AVE_inner)
        
    #a <- lapply(result$a, cbind)
    a = result$a
    d <- rep(0, J)
    for(j in seq_len(J)){
      d[j] <- correction(A[[j]], a[[j]], 
                         mode = mode[j], 
                         bias = bias)$d
    }
    
    names(d) = names(A)
    
    lambda = mapply("*", a, d, SIMPLIFY = FALSE)
    
    omega = list()
    for (j in which(mode == "formative")){
        omega[[j]] = solve(cov2(A[[j]], bias = bias))%*%a[[j]]*d[j]
    }
    
        
    for (b in 1:J) {
      rownames(a[[b]]) = colnames(A[[b]])  
      rownames(Y[[b]]) = rownames(A[[b]])
      colnames(Y[[b]]) = "comp1"
    }
        

    # reliability coefficients 
    reliability_coef = d^2/mapply("%*%", 
                                  lapply(a, t), 
                                  mapply("%*%", 
                                  lapply(A, function(x) cov2(x, bias= bias)), 
                                         a, SIMPLIFY = FALSE)
                                  )
        
      
    # rho_jk (Expected value : -1 <=rho_jk <=1)
    Phat = matrix(0, J, J)
    for (j in 1:J){
      for (k in j:J){
            Phat[j, k] = t(a[[j]])%*%cov2(A[[j]], A[[k]], bias = bias)%*%a[[k]]/(d[j]*d[k])
      }
    }
    
    Phat = Phat + t(Phat)
    diag(Phat) = 1
    colnames(Phat) = rownames(Phat) = names(A)
    beta_gamma = LV_MODEL(Phat,C)    
    
    # variance of the residuals
    # var_MVs = apply(Reduce("cbind", A), 2, 
    #                function(x) 
    #                cov2(x, bias = bias))
    
    var_MVs = lapply(A, function(x) diag(cov2(x, bias = bias)))
        
    residual_variance = mapply("-", var_MVs, lapply(lambda, function(x) x^2), 
                               SIMPLIFY = FALSE)
        
    # cor(y_jh, eta_j) 
    # Expected value : -1 <= cor(y_jh, eta_j) <=1
    cor_y_eta = lapply(1:length(lambda), 
                       function(x) lambda[[x]]/sqrt(var_MVs[[x]]))
    
    out <- list(a = a, 
                lambda = lambda,
                omega = omega,
                beta = beta_gamma[[1]],
                gamma = beta_gamma[[2]],
                d = d,
                crit = result$crit, 
                bias = bias,
                Phat = Phat,
                reliability_coef = reliability_coef,
                cor_y_eta = cor_y_eta,
                residual_variance = residual_variance, 
                C = C,
                mode = mode)
        
    class(out) <- "rgcca"
    return(out)
}
    
    