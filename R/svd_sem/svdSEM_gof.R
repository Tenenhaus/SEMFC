source("R/utils/scaleDataSet.R")

#global fit & test of fit
svdSEM_gof <- function(fit, B = 100, bias = FALSE){
  
  if(is.null(B)){
    return(T_LS = fit$T_LS)
  }else{
    # Tranforms the data sets in the way proposed 
    # by Yuan & Hayashi (2003)
    df <- Reduce("cbind", fit$blocks)
    Z0 <- scaleDataSet(df, fit$SIGMA_IMPLIED)
    Z0 <- lapply(split(data.frame(t(Z0)),
                      as.factor(rep(seq_along(fit$blocks),
                                    sapply(fit$blocks, NCOL)))),
                t)
    
    beta = fit$beta
    gamma = fit$gamma
    Tb_LS = rep(NA, B)
    
    Tb_LS <- pbapply::pbsapply(1:B,
                      function(b){
                        ind <- sample(NROW(df), replace = TRUE)
                        Zb <- lapply(Z0, function(x) x[ind, ])
                        S <- cov2(Reduce("cbind", Zb), bias = bias)
                        fit_b <- svdSEM(Zb,
                                       C = fit$C,
                                       scale = fit$scale, 
                                       mode = fit$mode, 
                                       bias = fit$bias)
                        if (!any(eigen(fit_b$Ptilde)$values<=0)){
                          return(NA)
                        }else{
                          return(d_LS(S, fit_b$SIGMA_IMPLIED))
                        }
                      }
                      )
    
    
    pval <- mean(fit$T_LS <= Tb_LS, na.rm = TRUE)
    return(list(T_LS = fit$T_LS,
                Tb_LS = Tb_LS,
                pval = pval, 
                improper = sum(is.na(Tb_LS))))
  }
}

  

