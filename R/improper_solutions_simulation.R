##################################################
############ Table IMPROPER SOLUTIONS ############
##################################################



N <- c(seq(20, 100, by = 10), 200, 300, 400, 500, 600, 700)
n_simu <-  1000
n_improper <- 9
n_improper_ml <- 7
improper_sol <- array(NA, dim = c(n_simu, n_improper, length(N)))
improper_sol_ml <- array(NA, dim = c(n_simu, n_improper_ml, length(N)))
Table_improper <- matrix(0, 9, length(N))
Table_improper_ml <- matrix(0, 7, length(N))


set.seed(27) #my favorite number

for (n in seq_along(N)){
  for(s in 1:n_simu){

    print(paste(N[n], " & ", s))
    X <- mvrnorm(N[n], rep(0, 18), SIGMA, empirical = FALSE)
    colnames(X) <- paste("X", rep(1:6, each = 3), rep(1:3, 6))

    Y <- list(X1 = X[, 1:3], X2 = X[, 4:6], X3 = X[, 7:9],
              X4 = X[, 10:12], X5 = X[, 13:15], X6 = X[, 16:18])


    model <- SemFC$new(data=Y, relation_matrix = C, mode=mode, scale=F, bias=F)
    model$fit_svd()


    improper_sol[s, , n] <- improper(model$parameters)

    dimnames(improper_sol) <- list(1:n_simu, c("RELIABILITY_COEF",
                                              "RHO_JH",
                                              "P_TILDE",
                                              "P_IMPLIED",
                                              "THETA_JH",
                                              "STD_LAMBDA",
                                              "SIGMA_IMPLIED",
                                              "R2",
                                              "PSI"), N)


    try({
      init_ml <- model$parameters$theta
      model_ml <- SemFC$new(data=Y, relation_matrix = C, mode=mode, scale=F, bias=F)
      model_ml$fit_ml(initialisation_svd = TRUE)

      # ML output

      x <- model_ml$parameters$theta

      SIGMA_ML <- model_ml$parameters$SIGMA_IMPLIED
      R_LVM_ML <- model_ml$parameters$R_LVM
      RESID_VAR_ML <- unlist(model_ml$parameters$residual_variance)
      STD_LAMBDA_ML <- unlist(model_ml$parameters$std_lambda)
      R2_ML <- unlist(model_ml$parameters$R2)
      PSI_ML <- model_ml$parameters$psi

      # Improper solution 1 : number of out of range rho_jk
      # Expected value : -1 <=rho_jk <=1

      improper_sol_ml[s, 1, n] <- any(abs(R_LVM_ML[upper.tri(R_LVM_ML)])>1)


      # Improper solution 2 : is P_IMPLIED positive definite
      # Expected value = TRUE

      improper_sol_ml[s, 2, n] <- any(eigen(R_LVM_ML)$values<0)


      # Improper solution 3 : is variance of the residual positive
      # Expected value = TRUE

      improper_sol_ml[s, 3, n] <- any(RESID_VAR_ML<0)

      # Improper solution 4 : number of out of range cor(y_jh, eta_j)
      # Expected value : -1 <= cor(y_jh, eta_j) <=1

      improper_sol_ml[s, 4, n] <- any(abs(STD_LAMBDA_ML)>1)

      # Improper solution 5 : is SIGMA_IMPLIED positive definite
      # Expected value = TRUE

      improper_sol_ml[s, 5, n] <- any(eigen(SIGMA_ML)$values<0)

      # Improper solution 6 : number of out of range R^2_i
      # Expected value : 0 <= R^2_i <=1

      improper_sol_ml[s, 6, n] <- any(R2_ML< 0 | R2_ML>1)

      # Improper solution 7 : is PSI positive definite
      # Expected value = TRUE

      improper_sol_ml[s, 7, n] <- any(eigen(PSI_ML)$values<0)

    }, silent = FALSE)


  }


}

dimnames(improper_sol_ml) <- list(1:n_simu, c("RHO_JH",
                                             "P_IMPLIED",
                                             "THETA_JH",
                                             "STD_LAMBDA",
                                             "SIGMA_IMPLIED",
                                             "R2",
                                             "PSI"), N)


Table_improper <- sapply(seq_along(N),
                         function(x) colMeans(improper_sol[, , x]))
colnames(Table_improper) <- N
Table_improper[-c(1, 3), ]


Table_improper_ml <- sapply(seq_along(N),
                            function(x) colMeans(improper_sol_ml[, , x]))
colnames(Table_improper_ml) <- N
Table_improper_ml
