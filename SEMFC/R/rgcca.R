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
#' @param C  A design matrix that describes the relationships between blocks (default: complete design).
#' @param tau tau is either a \eqn{1 * J} vector or a \eqn{max(ncomp) * J} matrix, and contains the values 
#' of the shrinkage parameters (default: tau = 1, for each block and each dimension).
#' If tau = "optimal" the shrinkage paramaters are estimated for each block and each dimension using the Schafer and Strimmer (2005)
#' analytical formula . If tau is a \eqn{1* J} numeric vector, tau[j] is identical across the dimensions of block \eqn{X_j}. 
#' If tau is a matrix, tau[k, j] is associated with \eqn{X_{jk}} (\eqn{k}th residual matrix for block \eqn{j})
#' @param ncomp  A \eqn{1 * J} vector that contains the numbers of components for each block (default: rep(1, length(A)), which gives one component per block.)
#' @param scheme The value is "horst", "factorial", "centroid" or any diffentiable convex scheme function g designed by the user (default: "factorial").
#' @param scale  If scale = TRUE, each block is standardized to zero means and unit variances and then divided by the square root of its number of variables (default: TRUE).
#' @param verbose  If verbose = TRUE, the progress will be report while computing (default: TRUE).
#' @param init The mode of initialization to use in RGCCA algorithm. The alternatives are either by Singular Value Decompostion ("svd") or random ("random") (Default: "svd").
#' @param bias A logical value for biaised or unbiaised estimator of the var/cov (default: bias = TRUE).
#' @param tol The stopping value for convergence.
#' @return \item{Y}{A list of \eqn{J} elements. Each element of \eqn{Y} is a matrix that contains the RGCCA components for the corresponding block.}
#' @return \item{a}{A list of \eqn{J} elements. Each element of \eqn{a} is a matrix that contains the outer weight vectors for each block.}
#' @return \item{astar}{A list of \eqn{J} elements. Each element of astar is a matrix defined as Y[[j]][, h] = A[[j]]\%*\%astar[[j]][, h].}
#' @return \item{C}{A design matrix that describes the relation between blocks (user specified).}
#' @return \item{tau}{A vector or matrix that contains the values of the shrinkage parameters applied to each block and each dimension (user specified).}
#' @return \item{scheme}{The scheme chosen by the user (user specified).}
#' @return \item{ncomp}{A \eqn{1 * J} vector that contains the numbers of components for each block (user specified).}
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


rgcca <- function(A, C = 1-diag(length(A)), tau = rep(1, length(A)), ncomp = rep(1, length(A)), 
                  scheme = "factorial", scale = TRUE , init="svd", bias = TRUE, tol = 1e-8, verbose=TRUE) {
    if (any(ncomp < 1)) stop("Compute at least one component per block!")
    pjs <- sapply(A, NCOL)
    nb_row <- NROW(A[[1]])
    if (any(ncomp-pjs > 0)) stop("For each block, choose a number of components smaller than the number of variables!")
    #-------------------------------------------------------
    if (mode(scheme) != "function") {
    if ((scheme != "horst" ) & (scheme != "factorial") & (scheme != "centroid")) {
      stop("Choose one of the three following schemes: horst, centroid, factorial or design the g function")
      }
    if (verbose) cat("Computation of the RGCCA block components based on the", scheme, "scheme \n")
    }
    if (mode(scheme) == "function" & verbose) {
      cat("Computation of the RGCCA block components based on the g scheme \n")
    }
    
    #-------------------------------------------------------
    
    if (scale == TRUE) {
      A = lapply(A, function(x) scale2(x, bias = bias))
      A = lapply(A, function(x) x/sqrt(NCOL(x)))
    }
    
    if (!is.numeric(tau) & verbose) {
      cat("Optimal Shrinkage intensity paramaters are estimated \n")
    } else {
      if (is.numeric(tau) & verbose) {
        cat("Shrinkage intensity paramaters are chosen manually \n")
      }
    }
    AVE_X = list()
    AVE_outer <- vector()
    lambda = list()
    ndefl <- ncomp-1
    N <- max(ndefl)
    nb_ind <- NROW(A[[1]])
    J <- length(A)
    primal_dual = rep("primal", J)
    primal_dual[which(nb_row<pjs)] = "dual"
    improper_solution = list()
    
    if (N == 0) {
        result <- rgccak(A, C, tau=tau , scheme=scheme, init = init, bias = bias, tol = tol, verbose=verbose)
        Y <- NULL
        for (b in 1:J) Y[[b]] <- result$Y[,b, drop = FALSE]
        
        #Average Variance Explained (AVE) per block
        for (j in 1:J) AVE_X[[j]] =  mean(cor(A[[j]], Y[[j]])^2)
        
        #AVE outer 
        AVE_outer <- sum(pjs * unlist(AVE_X))/sum(pjs)
                
        AVE <- list(AVE_X = AVE_X, 
                    AVE_outer = AVE_outer,
                    AVE_inner = result$AVE_inner)
        
        a <- lapply(result$a, cbind)
        
        disattenuation = mapply(function(x, y) correction(x, y, bias = bias), 
                                A, a, SIMPLIFY = TRUE)
        
        d = unlist(disattenuation[1, ])
        
        lambda = mapply("*", a, d, SIMPLIFY = FALSE)
        
        for (b in 1:J) {
          rownames(a[[b]]) = colnames(A[[b]])  
          rownames(Y[[b]]) = rownames(A[[b]])
          colnames(Y[[b]]) = "comp1"
        }
        
        ######################
        # IMPROPER SLOUTIONS #
        ######################
        
        # Improper solution 1 : number of negative d_jh^2
        improper_solution[[1]] = unlist(disattenuation[2, ])
        
        # Improper solution 2 : number of out of range reliability coefficients 
        # Expected value : 0 <=rho_jj <=1
        reliability_coef = d^2/mapply("%*%", 
                                    lapply(a, t), 
                                    mapply("%*%", 
                                           lapply(A, function(x) cov2(x, bias= bias)), 
                                           a, SIMPLIFY = FALSE)
                                    )
        
        range_reliability_coef = rep(TRUE, J)
        range_reliability_coef[which(reliability_coef>0 & reliability_coef<1)] = FALSE
        improper_solution[[2]] = range_reliability_coef
        
        # Improper solution 3 : number of out of range rho_jk 
        # Expected value : -1 <=rho_jk <=1
        
        Phat = matrix(0, J, J)
        for (j in 1:J){
          for (k in j:J){
            Phat[j, k] = t(a[[j]])%*%cov2(A[[j]], A[[k]], bias = bias)%*%a[[k]]/(d[j]*d[k])
          }
        }
        Phat = Phat + t(Phat)
        diag(Phat) = 1
        
        improper_solution[[3]] = sum(abs(Phat[upper.tri(Phat)])>1)
        
        # Improper solution 4 : is Phat non positive definite 
        # Expected value = TRUE
        
        improper_solution[[4]] = any(eigen(Phat)$values<0)

        # Improper solution 5 : is variance of the residual positive
        # Expected value = TRUE
        
        var_MVs = apply(Reduce("cbind", A), 2, 
                        function(x) 
                          cov2(x, bias = bias))
        
        residual_variance = var_MVs - unlist(lambda)^2
        
        improper_solution[[5]] = sum(residual_variance<0)
        
        # Improper solution 6 : number of out of range cor(y_jh, eta_j) 
        # Expected value : -1 <= cor(y_jh, eta_j) <=1
        
        cor_y_eta = unlist(lambda)/sqrt(var_MVs)
        
        improper_solution[[6]] = sum(abs(cor_y_eta)>1)
        
        names(improper_solution) = c("d_jk<0", "oor_reliability_coef",
                                     "oor_rhojk", "Phat<0", 
                                     "theta_jh<0", "oor_cor_yjh_eta")
        
        
        
        out <- list(Y=Y, 
                    a =a, 
                    lambda = lambda, 
                    astar=a, 
                    AVE=AVE, 
                    d = d,
                    crit=result$crit, 
                    C=C, tau=result$tau, scheme=scheme, 
                    ncomp=ncomp, primal_dual = primal_dual,
                    reliability_coef = reliability_coef,
                    cor_y_eta = cor_y_eta,
                    Phat = Phat,
                    improper_solution = improper_solution
                    )
        
        class(out) <- "rgcca"
        return(out)
    }
    
    Y <- NULL
    crit = d = list()
    AVE_inner <- rep(NA,max(ncomp))
    
    R <- A
    P <- a <- astar <- lambda <- NULL
    d = matrix(0, max(ncomp), J)
    neg_dk = matrix(0, max(ncomp), J)
    
    if (is.numeric(tau)) tau_mat = tau
    else tau_mat = matrix(NA, max(ncomp), J)
    
    for (b in 1:J) P[[b]] <- a[[b]] <- astar[[b]] <- lambda[[b]] <- matrix(NA,pjs[[b]],N+1)
    for (b in 1:J) Y[[b]] <- matrix(NA,nb_ind,N+1)
    
    for (n in 1:N) {
        if (verbose) cat(paste0("Computation of the RGCCA block components #", n, " is under progress...\n"))
        if(is.vector(tau)) 
          rgcca.result <- rgccak(R, C, tau = tau , scheme=scheme, init = init, bias = bias, tol = tol, verbose=verbose) 
        else
          rgcca.result <- rgccak(R, C, tau = tau[n, ] , scheme=scheme, init = init, bias = bias, tol = tol, verbose=verbose) 
        
        if (!is.numeric(tau)) tau_mat[n, ] = rgcca.result$tau
        AVE_inner[n] <- rgcca.result$AVE_inner
        crit[[n]] <- rgcca.result$crit
        
        d[n, ] = unlist(mapply(function(x, y) correction(x, y, bias = bias), 
                               R, rgcca.result$a)[1, ])
        
        neg_djk[n, ] = unlist(mapply(function(x, y) correction(x, y, bias = bias), 
                                     R, rgcca.result$a)[2, ])
        
        for (b in 1:J) {
          lambda[[b]][,n] <- rgcca.result$a[[b]]*d[n, b]
        }
        
        for (b in 1:J) Y[[b]][,n] <- rgcca.result$Y[ , b]
	      defla.result <- defl.select(rgcca.result$Y, R, ndefl, n, nbloc = J)
        R <- defla.result$resdefl
        for (b in 1:J) P[[b]][,n] <- defla.result$pdefl[[b]]
        for (b in 1:J) a[[b]][,n] <- rgcca.result$a[[b]]
        
        
        if (n==1) {
           for (b in 1:J) astar[[b]][,n] <- rgcca.result$a[[b]]
        } 
        else {
           for (b in 1:J) astar[[b]][,n] <- rgcca.result$a[[b]] - astar[[b]][,(1:n-1), drop=F] %*% drop( t(a[[b]][,n]) %*% P[[b]][,1:(n-1),drop=F] ) 
        }
    }
    
    if (verbose) cat(paste0("Computation of the RGCCA block components #", N+1, " is under progress ... \n"))
    if(is.vector(tau)) 
      rgcca.result <- rgccak(R, C, tau = tau , scheme=scheme, init = init, bias = bias, tol = tol, verbose=verbose) 
    else
      rgcca.result <- rgccak(R, C, tau = tau[N+1, ] , scheme=scheme, init = init, bias = bias, tol = tol, verbose=verbose) 
    
    crit[[N+1]] <- rgcca.result$crit

    if (!is.numeric(tau)) tau_mat[N+1, ] = rgcca.result$tau    
    AVE_inner[max(ncomp)] <- rgcca.result$AVE_inner
    
    d[N+1, ] = mapply(function(x, y) correction(x, y, bias = bias), R, rgcca.result$a)
    
    for (b in 1:J) {
      lambda[[b]][,N+1] <- rgcca.result$a[[b]]*d[N+1, b]
    }
    
    for (b in 1:J) {
      Y[[b]][,N+1]     <- rgcca.result$Y[ ,b]
      a[[b]][,N+1]     <- rgcca.result$a[[b]]
      astar[[b]][,N+1] <- rgcca.result$a[[b]] - astar[[b]][,(1:N),drop=F] %*% drop( t(a[[b]][,(N+1)]) %*% P[[b]][,1:(N),drop=F] ) 
      rownames(a[[b]]) = rownames(astar[[b]]) = rownames(lambda[[b]]) = colnames(A[[b]])  
      rownames(Y[[b]]) = rownames(A[[b]])
      colnames(Y[[b]]) = paste0("comp", 1:max(ncomp))
    }
    
    colnames(d) = paste0("Block", 1:J)
    rownames(d) = paste0("comp", 1:max(ncomp))
    
    shave.matlist <- function(mat_list, nb_cols) mapply(function(m, nbcomp) m[, 1:nbcomp, drop = FALSE], mat_list, nb_cols,SIMPLIFY=FALSE)
    shave.veclist <- function(vec_list, nb_elts) mapply(function(m, nbcomp) m[1:nbcomp], vec_list, nb_elts, SIMPLIFY=FALSE)    
    
    #Average Variance Explained (AVE) per block
    for (j in 1:J) AVE_X[[j]] =  apply(cor(A[[j]], Y[[j]])^2, 2, mean)
    
    #AVE outer 
    outer = matrix(unlist(AVE_X), nrow = max(ncomp))
    for (j in 1:max(ncomp)) AVE_outer[j] <- sum(pjs * outer[j, ])/sum(pjs)
    
    Y = shave.matlist(Y, ncomp)
    AVE_X = shave.veclist(AVE_X, ncomp)
    
    AVE <- list(AVE_X = AVE_X, 
                AVE_outer_model = AVE_outer,
                AVE_inner_model = AVE_inner)
    
    out <- list(Y = shave.matlist(Y, ncomp),
                a = shave.matlist(a, ncomp), 
                lambda = lambda,
                astar = shave.matlist(astar, ncomp),
                C = C, tau = tau_mat, 
                scheme = scheme, ncomp=ncomp, crit = crit,
                primal_dual = primal_dual,
                AVE = AVE, 
                d = d,
                improper_solution = improper_solution,
                )
    
    class(out) <- "rgcca"
    return(out)
}
