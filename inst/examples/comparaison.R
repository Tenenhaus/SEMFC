



F1_oldvers <- function(x, S){
  #loadings
  l1 = x[1:3] ; l2 = x[4:6] ; l3 = x[7:9]
  l4 = x[10:12] ; l5 = x[13:15] ;   l6 = x[16:18]
  L = bdiag(list(l1, l2, l3, l4, l5, l6))
  # correlations between exogeneous
  P_EXO = matrix(c(1, x[19], x[20], x[22],
                   x[19], 1, x[21], x[23],
                   x[20], x[21], 1, x[24],
                   x[22], x[23], x[24], 1), 4, 4)
  # path coefficients
  G = matrix(c(x[25], x[26], 0, 0,
               0, 0, x[27], x[28]), 2, 4, byrow = TRUE)

  # path coefficients
  B = matrix(c(0, x[29],
               x[30], 0), 2, 2, byrow = TRUE)

  # Correlations between endogeneous
  P_ENDO = matrix(c(1, x[31],
                    x[31], 1), 2, 2)

  # Correlations between Latent/emergent variables
  R = rbind(cbind(P_EXO, P_EXO%*%t(G)%*%t(solve(diag(2) - B))),
            cbind(solve(diag(2) - B)%*%G%*%P_EXO, P_ENDO))

  #cov between MVs
  S1 = matrix(c(x[32], x[33], x[35],
                x[33], x[34], x[36],
                x[35], x[36], x[37]), 3, 3)
  #cov between MVs
  S2 = matrix(c(x[38], x[39], x[41],
                x[39], x[40], x[42],
                x[41], x[42], x[43]), 3, 3)
  #cov between MVs
  S3 = matrix(c(x[44], x[45], x[47],
                x[45], x[46], x[48],
                x[47], x[48], x[49]), 3, 3)
  #cov between MVs
  S4 = matrix(c(x[50], x[51], x[53],
                x[51], x[52], x[54],
                x[53], x[54], x[55]), 3, 3)

  T5 = diag(x[56:58])
  T6 = diag(x[59:61])

  implied_S = L%*%R%*%t(L) +
    bdiag(list(S1 - l1%*%t(l1),
               S2 - l2%*%t(l2),
               S3 - l3%*%t(l3),
               S4 - l4%*%t(l4),
               T5,
               T6))

  opt = log(det(implied_S)) +
    sum(diag(S%*%solve(implied_S))) -
    log(det(S)) -
    NCOL(S)

  return(opt)
}


heq1_oldvers <- function(x, S) {

  #cov between MVs
  S1 = matrix(c(x[32], x[33], x[35],
                x[33], x[34], x[36],
                x[35], x[36], x[37]), 3, 3)
  #cov between MVs
  S2 = matrix(c(x[38], x[39], x[41],
                x[39], x[40], x[42],
                x[41], x[42], x[43]), 3, 3)
  #cov between MVs
  S3 = matrix(c(x[44], x[45], x[47],
                x[45], x[46], x[48],
                x[47], x[48], x[49]), 3, 3)
  #cov between MVs
  S4 = matrix(c(x[50], x[51], x[53],
                x[51], x[52], x[54],
                x[53], x[54], x[55]), 3, 3)

  l1 = x[1:3] ; l2 = x[4:6] ; l3 = x[7:9]
  l4 = x[10:12]

  h <- c(rep(0,4))
  h[1] <- t(l1)%*%solve(S1)%*%l1 - 1
  h[2] <- t(l2)%*%solve(S2)%*%l2 - 1
  h[3] <- t(l3)%*%solve(S3)%*%l3 - 1
  h[4] <- t(l4)%*%solve(S4)%*%l4 - 1
  return(h)
}



set.seed(1)
source('R/SEMFC/sem_f_c.R')

source('inst/model/model_mixed.R')

model <- SemFC$new(data=Y, relation_matrix = C, mode=mode, scale=F, bias=F)
model$fit_svd()




init_ml_with_S <- model$parameters$theta


####################################  NEW VERSION  ####################################

r <- sum(mode == "formative")



F1_new = F1(x = init_ml_with_S,
            S = model$cov_S, block_sizes=model$block_sizes, mode=model$mode, lengths_parameter = model$lengths_theta,
            which_exo_endo = model$which_exo_endo)


implied_S_new <- lvm_ml(init_ml_with_S, model$block_sizes, model$mode, model$lengths_theta, model$which_exo_endo, jac = FALSE)

heq1_new = heq1(x = init_ml_with_S,
                  S = model$cov_S, block_sizes=model$block_sizes, mode=model$mode, lengths_parameter = model$lengths_theta,
                  which_exo_endo = model$which_exo_endo)



result_new <- solnp(pars = init_ml_with_S,
                fun=F1, eqfun=heq1,
                eqB = rep(0,r), S = model$cov_S, block_sizes=model$block_sizes, mode=model$mode, lengths_parameter = model$lengths_theta,
                which_exo_endo = model$which_exo_endo,
                control = list(trace = 0, tol = 1e-8))


x.new = result_new$pars

model_ml <- SemFC$new(data=Y, relation_matrix = C, mode=mode, scale=F, bias=F)
model_ml$fit_ml(initialisation_svd  = TRUE)

x.ml <- model_ml$parameters$theta

implied_S <- lvm_ml(init_ml_with_S, model$block_sizes, model$mode, model$lengths_theta, model$which_exo_endo, jac = FALSE)$SIGMA_IMPLIED

####################################  OLD VERSION  ####################################





F1_old = F1_oldvers(x = init_ml_with_S,
                      S = model$cov_S)

heq1_old = heq1_oldvers(x = init_ml_with_S,
                            S = model$cov_S)




result_old = solnp(pars = init_ml_with_S,
                     fun=F1_oldvers, eqfun=heq1_oldvers,
                     eqB = rep(0,4), S = model$cov_S,
                     control = list(trace = 0, tol = 1e-8))


x.old = result_old$pars

x  =init_ml_with_S

P_EXO = matrix(c(1, x[19], x[20], x[22],
                 x[19], 1, x[21], x[23],
                 x[20], x[21], 1, x[24],
                 x[22], x[23], x[24], 1), 4, 4)

P_ENDO = matrix(c(1, x[31],
                  x[31], 1), 2, 2)

# path coefficients
G = matrix(c(x[25], x[26], 0, 0,
             0, 0, x[27], x[28]), 2, 4, byrow = TRUE)

# path coefficients
B = matrix(c(0, x[29],
             x[30], 0), 2, 2, byrow = TRUE)

R_LVM_ML = rbind(cbind(P_EXO, P_EXO%*%t(G)%*%t(solve(diag(2) - B))),
                 cbind(solve(diag(2) - B)%*%G%*%P_EXO, P_ENDO))

I_B_ML = diag(2) - B

PSI_ML =  I_B_ML%*%P_ENDO%*%t(I_B_ML) - G%*%P_EXO%*%t(G)

r2_1_ML = 1 - PSI_ML[1, 1]
r2_2_ML = 1 - PSI_ML[2, 2 ]

#Omega
#cov between MVs
S1 = matrix(c(x[32], x[33], x[35],
              x[33], x[34], x[36],
              x[35], x[36], x[37]), 3, 3)

#cov between MVs

S2 = matrix(c(x[38], x[39], x[41],
              x[39], x[40], x[42],
              x[41], x[42], x[43]), 3, 3)

#cov between MVs
S3 = matrix(c(x[44], x[45], x[47],
              x[45], x[46], x[48],
              x[47], x[48], x[49]), 3, 3)
#cov between MVs
S4 = matrix(c(x[50], x[51], x[53],
              x[51], x[52], x[54],
              x[53], x[54], x[55]), 3, 3)

LAMBDA_ML = bdiag(split(x[1:18], f = rep(1:6, each = 3)))
RESID_VAR_ML = apply(X, 2, var) - x[1:18]^2
SIGMA_ML = LAMBDA_ML%*%R_LVM_ML%*%t(LAMBDA_ML) + diag(RESID_VAR_ML)

SIGMA_ML[1:3, 1:3] = S1      #S[1:3, 1:3]
SIGMA_ML[4:6, 4:6] = S2      #S[4:6, 4:6]
SIGMA_ML[7:9, 7:9] = S3      #S[7:9, 7:9]
SIGMA_ML[10:12, 10:12] = S4  #S[10:12, 10:12]






