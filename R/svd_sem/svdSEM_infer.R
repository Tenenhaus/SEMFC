svdSEM_infer <- function(fit, B = 100, verbose = TRUE){
  
  df <- data.frame(Reduce("cbind", fit$blocks))
  sd_init <- apply(df, 2, sd)
  
  beta <- fit$beta[fit$beta!=0]
  gamma <- fit$gamma[fit$gamma!=0]
  lambda <- unlist(unname(fit$lambda))
  # names(lambda) <- unlist(lapply(fit$lambda, names))
  std_loadings <- lambda/sd_init
  residual_variance <- unlist(unname(fit$residual_variance))
  
  Z0 <- lapply(split(data.frame(t(df)),
                    as.factor(rep(seq_along(fit$blocks), 
                                  sapply(fit$blocks, NCOL)))),
              t)
  
  if(verbose == FALSE){
    pbapply::pboptions(type = "none") 
  }else{
    pbapply::pboptions(type = "txt")
  }
  
  L <- pbapply::pbsapply(1:B,
                        function(b){
                          ind <- sample(NROW(df), replace = TRUE)
                          Zb <- lapply(Z0, function(x) x[ind, ])
                          names(Zb) <- names(fit$blocks)
                          SD <- apply(Reduce("cbind", Zb), 2, sd)
                          fit_b <- svdSEM(Zb,
                                         C = fit$C, 
                                         scale = fit$scale, 
                                         mode = fit$mode, 
                                         bias = fit$bias)
                          if (!any(eigen(fit_b$Ptilde)$values<=0)){
                            return(list(
                              boot_lambda  = as.vector(Reduce("c", fit_b$lambda)),
                              boot_std_loadings = as.vector(Reduce("c", fit_b$lambda)/SD),
                              boot_beta = fit_b$beta[fit_b$beta!=0],
                              boot_gamma = fit_b$gamma[fit_b$gamma!=0],
                              boot_residual_variance = as.vector(Reduce("c", fit_b$residual_variance)))
                            )
                          }
                          else{
                            return(list(
                            boot_lambda = NA,
                            boot_std_loadings = NA,
                            boot_beta = NA,
                            boot_gamma = NA,
                            boot_residual_variance = NA)
                            )
                          }
                        } 
                        )
  
  boot_lambda <- Reduce("rbind", L[1, ][!is.na(L[1, ])])
  boot_std_loadings <- Reduce("rbind", L[2, ][!is.na(L[2, ])])
  boot_beta <- Reduce("rbind", L[3, ][!is.na(L[3, ])])
  boot_gamma <- Reduce("rbind", L[4, ][!is.na(L[4, ])])
  boot_residual_variance <- Reduce("rbind", L[5, ][!is.na(L[5, ])])


  std_residual_variance <- apply(boot_residual_variance, 2, sd)
  t_ratio <- residual_variance/std_residual_variance
  pval_residual_variance <- sapply(1:length(t_ratio),
                             function(x)
                               2*pnorm(abs(t_ratio[x]),
                                       lower.tail = FALSE)
  )
  residual_variance <- data.frame(lambda = residual_variance,
                      std = std_residual_variance,
                      z = t_ratio,
                      pval = pval_residual_variance)


  
  std_lambda <- apply(boot_lambda, 2, sd)
  
  t_ratio <- lambda/std_lambda
  pval_lambda <- sapply(1:length(t_ratio),
                             function(x)
                               2*pnorm(abs(t_ratio[x]),
                                       lower.tail = FALSE)
  )
  
  
  lambda <- data.frame(lambda = lambda,
                      std = std_lambda,
                      z = t_ratio,
                      pval = pval_lambda)
  
  std_std_loadings <- apply(boot_std_loadings, 2, sd)
  
  t_ratio <- std_loadings/std_std_loadings
  pval_std_loadings <- sapply(1:length(t_ratio),
                              function(x)
                                2*pnorm(t_ratio[x],
                                        lower.tail = FALSE)
  )
  
  
  std_lambda = data.frame(std_loadings = std_loadings, 
                          std = std_std_loadings,
                          z = t_ratio,
                          pval = pval_std_loadings)
  
  std_beta = apply(boot_beta, 2, sd)
  t_ratio = beta[beta!=0]/apply(boot_beta, 2, sd)
  pval_beta = sapply(1:length(t_ratio),
                         function(x)
                           2*pnorm(abs(t_ratio[x]), lower.tail = FALSE)
  )

  beta <- data.frame()
  if (length(beta) != 0){

    beta = data.frame(beta = beta,
                    std = std_beta,
                    z = t_ratio,
                    pval = pval_beta)

    rownames(beta) = sapply(1:NROW(beta),
           function(b)
             paste(colnames(fit$beta)[which(fit$beta!=0, arr.ind = TRUE)[b, ]],
                   collapse = "~")
           )

  }


  std_gamma = apply(boot_gamma, 2, sd)
  
  t_ratio = gamma[gamma!=0]/std_gamma
  pval_gamma = sapply(1:length(t_ratio),
                          function(x)
                            2*pnorm(abs(t_ratio[x]), lower.tail = FALSE)
  )
  
  gamma = data.frame(gamma = gamma, 
                     std = std_gamma,
                     z = t_ratio,
                     pval = pval_gamma)
  
  rownames(gamma) = sapply(1:NROW(gamma),
                          function(b) 
                            paste(rownames(fit$gamma)[which(fit$gamma!=0, arr.ind = TRUE)[b, 1]], 
                                  colnames(fit$gamma)[which(fit$gamma!=0, arr.ind = TRUE)[b, 2]],
                                  sep = "~")
  )
  
  
  return(list(out = list(boot_lambda, 
                         boot_std_loadings,
                         boot_beta,
                         boot_gamma,
                         boot_residual_variance),
              lambda = lambda,
              std_lambda = std_lambda,
              beta = beta,
              gamma = gamma,
              residual_variance = residual_variance,
              improper = sum(is.na(L[1, ]))))
}
  

