improper <- function(fit){
  
  improper_sol = rep(FALSE, 9)
  
  # Improper solution 1 : number of out of range reliability coefficient
  # Expected value : 0 <=rho_jk <=1

  improper_sol[1] = any(fit.svd$reliability_coef < 0 | 
                        fit.svd$reliability_coef > 1)
  
  
  # Improper solution 2 : number of out of range rho_jk 
  # Expected value : -1 <=rho_jk <=1
  
  improper_sol[2] = any(abs(fit.svd$P_IMPLIED[upper.tri(fit.svd$P_IMPLIED)])>1)
  
  # Improper solution 3 : is P_TILDE positive definite 
  # Expected value = TRUE
  
  improper_sol[3] = any(eigen(fit.svd$Ptilde)$values<0)
  
  # Improper solution 4 : is P_IMPLIED positive definite 
  # Expected value = TRUE
  
  improper_sol[4] = any(eigen(fit.svd$P_IMPLIED)$values<0)
  
  # Improper solution 5 : is variance of the residual positive
  # Expected value = TRUE
  
  improper_sol[5] = any(unlist(fit.svd$residual_variance)<0)
  
  # Improper solution 6 : number of out of range cor(y_jh, eta_j) 
  # Expected value : -1 <= cor(y_jh, eta_j) <=1
  
  improper_sol[6] = any(abs(unlist(fit.svd$std_lambda))>1)
  
  # Improper solution 7 : is SIGMA_IMPLIED positive definite 
  # Expected value = TRUE
  
  improper_sol[7] = any(eigen(fit.svd$SIGMA_IMPLIED)$values<0)
  
  # Improper solution 8 : number of out of range R^2_i
  # Expected value : 0 <= R^2_i <=1
  
  improper_sol[8] = any(fit.svd$R2< 0 | fit.svd$R2>1)
  
  # Improper solution 9 : is PSI positive definite 
  # Expected value = TRUE
  
  improper_sol[9] = any(eigen(fit.svd$psi)$values<0) 
  
  names(improper_sol) = c("RELIABILITY_COEF",
                          "RHO_JH", 
                          "P_TILDE", 
                          "P_IMPLIED", 
                          "THETA_JH", 
                          "STD_LAMBDA", 
                          "SIGMA_IMPLIED",
                          "R2",
                          "PSI")
  
  return(improper_sol= improper_sol)
}