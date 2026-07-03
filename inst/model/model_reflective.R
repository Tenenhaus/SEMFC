source('inst/simulations/data_simulation_reflective.R')

library(MASS)
library(mvtnorm)

set.seed(20091979)
N <- 300

X <- mvrnorm(N, rep(0, 18), SIGMA, empirical = TRUE)
colnames(X) <- paste("X", rep(1:6, each = 3), rep(1:3, 6), sep ="")


Y <- list(LV1 = X[, 1:3], LV2 = X[, 4:6], LV3 = X[, 7:9],
         LV4 = X[, 10:12], LV5 = X[, 13:15], LV6 = X[, 16:18])



C <- matrix(c(0, 0, 0, 0, 1, 0,
             0, 0, 0, 0, 1, 0,
             0, 0, 0, 0, 0, 1,
             0, 0, 0, 0, 0, 1,
             0, 0, 0, 0, 0, 1,
             0, 0, 0, 0, 1, 0), 6, 6, byrow = TRUE)

colnames(C) <- rownames(C) <- names(Y)

mode <- rep("reflective", 6)



X_2 <- mvrnorm(N, rep(0, 18), SIGMA, empirical = FALSE)
colnames(X_2) <- paste("X", rep(1:6, each = 3), rep(1:3, 6), sep ="")


Y_2 <- list(LV1 = X_2[, 1:3], LV2 = X_2[, 4:6], LV3 = X_2[, 7:9],
         LV4 = X_2[, 10:12], LV5 = X_2[, 13:15], LV6 = X_2[, 16:18])





sem.model <-  '
# latent variable definitions
LV1 =~ X11+X12+X13
LV2 =~ X21+X22+X23
LV3 =~ X31+X32+X33
LV4 =~ X41+X42+X43
LV5 =~ X51+X52+X53
LV6 =~ X61+X62+X63

# Regressions
LV5 ~ LV1 + LV2 + LV6
LV6 ~ LV3 + LV4 + LV5

# residual covariances
LV5 ~~ LV6
'






