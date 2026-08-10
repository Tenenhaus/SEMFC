source('inst/simulations/data_simulation_mixed.R')

library(MASS)
library(mvtnorm)
# set.seed(20091979)


X <- mvrnorm(N, rep(0, 6*q), SIGMA, empirical = TRUE)
colnames(X) <- paste("X", rep(1:6, each = q), rep(1:q, 6), sep ="")


Y <- list(LV1 = X[, 1:3], LV2 = X[, 4:6], LV3 = X[, 7:9],
         LV4 = X[, 10:12], LV5 = X[, 13:15], LV6 = X[, 16:18])



C <- matrix(c(0, 0, 0, 0, 1, 0,
             0, 0, 0, 0, 1, 0,
             0, 0, 0, 0, 0, 1,
             0, 0, 0, 0, 0, 1,
             0, 0, 0, 0, 0, 1,
             0, 0, 0, 0, 1, 0), 6, 6, byrow = TRUE)

colnames(C) <- rownames(C) <- names(Y)

mode <- c(rep("formative", 4), rep("reflective", 2))


X_2 <- mvrnorm(N, rep(0, q*6), SIGMA, empirical = FALSE)
colnames(X_2) <- paste("X", rep(1:6, each = q), rep(1:q, 6), sep ="")


Y_2 <- list(LV1 = X_2[, 1:q], LV2 = X_2[, (q+1):(2*q)], LV3 = X_2[, (2*q+1):(3*q)],
         LV4 = X_2[, (3*q+1):(4*q)], LV5 = X_2[, (4*q+1):(5*q)], LV6 = X_2[, (5*q+1):(6*q)])



sem.model <-  '
# latent variable definitions
eta5 =~ X51+X52+X53
eta6 =~ X61+X62+X63

# Composite model
eta1 <~ X11 + X12 + X13
eta2 <~ X21 + X22 + X23
eta3 <~ X31 + X32 + X33
eta4 <~ X41 + X42 + X43

# Regressions
eta5 ~ eta1 + eta2 + eta6
eta6 ~ eta3 + eta4 + eta5

# residual covariances
eta5 ~~ eta6
'


