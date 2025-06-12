diagnostic <- function(A, obj.rgccac){
  
  improper_solution = list()
  
######################
# IMPROPER SOLUTIONS #
######################
a <- lapply(obj.rgccac$a, cbind)
d <- obj.rgccac$d
bias <- obj.rgccac$bias
lambda <- obj.rgccac$lambda  
    
  ##################################################
  # Improper solution 1: number of negative d_jk^2 #
  ##################################################
  
    disattenuation = mapply(function(x, y) 
                              correction(x, y, bias = bias),
                            A, a, SIMPLIFY = TRUE)
  
    improper_solution[[1]] <- as.numeric(any(unlist(disattenuation[2, ])!=0)) 
      
  
  ###############################################################
  # Improper solution 2 : out of range reliability coefficients #
  #              Expected value : 0 <=rho_jj <=1                #
  ###############################################################
    
    reliability_coef <- obj.rgccac$reliability_coef
    range_reliability_coef = rep(TRUE, J)
    range_reliability_coef[which(reliability_coef>0 & reliability_coef<1)] = FALSE
    improper_solution[[2]] = as.numeric(any(range_reliability_coef!=0))
    
  #######################################################
  # Improper solution 3 : number of out of range rho_jk #
  #          Expected value : -1 <=rho_jk <=1           #
  #######################################################
    
    Phat <- fit.rgccac$Phat
    improper_solution[[3]] = as.numeric(any(abs(Phat[upper.tri(Phat)])>1))
  
  #######################################################
  # Improper solution 4 : is Phat non positive definite #
  #               Expected value = TRUE                 #
  #######################################################
  
    improper_solution[[4]] = as.numeric(any(eigen(Phat)$values<0))
 
  ##############################################################
  # Improper solution 5 : is variance of the residual positive #
  #               Expected value = TRUE                        #
  ##############################################################
  
    var_MVs = apply(Reduce("cbind", A), 2, 
                    function(x) 
                      cov2(x, bias = bias)
                    )
  
    residual_variance = var_MVs - unlist(lambda)^2
    improper_solution[[5]] = as.numeric(any(residual_variance<0))

    
  #################################################################
  # Improper solution 6 : number of out of range cor(y_jh, eta_j) #
  #       Expected value : -1 <= cor(y_jh, eta_j) <=1             #
  #################################################################
  
  cor_y_eta = unlist(lambda)/sqrt(var_MVs)
  improper_solution[[6]] = as.numeric(any(abs(cor_y_eta)>1))
  
  names(improper_solution) = c("d_jk", "reliability_coef",
                               "rho_jk", "Phat", 
                               "theta_jh", "cor_yjh_eta")
  
  
 return(improper_solution) 
  
}