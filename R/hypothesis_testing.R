

#########################################
###      hypothesis testing with      ###
###           ML-SEM, eigenSEM        ###
#########################################


# Bootstrap testing
nsimu <- 1000
nboot <- 1000
N <- c(300, 600, 1200)



decision_svd <- matrix(NA, nsimu, length(N))
decision_ml <- matrix(NA, nsimu, length(N))

T_LS_SVD <- matrix(NA, nsimu, length(N))
Tb_LS_SVD <- array(NA, c(nboot, nsimu, length(N)))
f_ml <- matrix(NA, nsimu, length(N))

set.seed(20091979)

for(n in seq_along(N)){

  for (s in 1:nsimu){
    X <- mvrnorm(N[n], rep(0, 18), SIGMA, empirical = FALSE)
    colnames(X) <- paste0("X", rep(1:6, each = 3), rep(1:3, 6))


    Y <- list(X1 = X[, 1:3], X2 = X[, 4:6], X3 = X[, 7:9],
         X4 = X[, 10:12], X5 = X[, 13:15], X6 = X[, 16:18])

    ##########
    #  SVD   #
    ##########


    model <- SemFC$new(data=Y, relation_matrix = C, mode=mode, scale=F, bias=F)
    model$fit_svd()

    P_LVM_SVD <- model$parameters$P_IMPLIED
    eigenval <- eigen(P_LVM_SVD)$values


    if(!any(eigenval<0)){
      T_LS_SVD[s, n] <- model$parameters$T_LS

      try({

        ##########
        #   ML   #
        ##########
        model_ml <- SemFC$new(data=Y, relation_matrix = C, mode=mode, scale=F, bias=F)
        model_ml$fit_ml(initialisation_svd = TRUE)


        f_ml[s, n] <- model_ml$parameters$F


      })


    }
    # transform the data sets in the way proposed
    # by Yuan & Hayashi (2003)

    if(!any(eigenval<0)){
      Z_SVD <- scaleDataSet(X, model$parameters$SIGMA_IMPLIED)
      Z_SVD <- lapply(split(data.frame(t(Z_SVD)),
                           as.factor(rep(1:6, each = 3))), t)
    }

    ################################
    # Bootstrap hypotestis testing #
    ################################

    for (b in 1:nboot){
      ind <- sample(1:N[n], replace = TRUE)

      ##########
      # SVD  #
      ##########

      if(!any(eigenval<0)){
        Zb_SVD <- lapply(Z_SVD, function(x) x[ind, ])
        Sb_SVD <- cov2(Reduce("cbind", Zb_SVD), bias = FALSE)

        model <- SemFC$new(data=Y, relation_matrix = C, mode=mode, scale=F, bias=F)
        model$fit_svd()

        # SIGMA_SEM_SVD = fit.svd$SIGMA_IMPLIED
        Tb_LS_SVD[b, s, n] <- model$parameters$T_LS
      }

    }
    decision_ml[s, n] <- pchisq((N[n]-1)*f_ml[s, n], 114, lower.tail = FALSE)<=0.05
    decision_svd[s, n] <- T_LS_SVD[s, n] > quantile(Tb_LS_SVD[, s, n], 0.95, na.rm = TRUE)

  }

}


######################
####### Table 7 ######
# global test-of-fit #
######################
######################



Table7 <- matrix(NA, 2, 3)

Table7 <- rbind(colMeans(decision_svd),
                    colMeans(decision_ml))

rownames(Table7) <- c("SVD-SEM", "ML-SEM")
colnames(Table7) <- c("N=300", "N=600", "N=1200")
round(Table7, 3)





