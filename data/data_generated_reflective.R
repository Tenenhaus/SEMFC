source('data/data_simulation_reflective.R')


# set.seed(20091979)
N <- 300
library(mvtnorm)
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
eta1 =~ X11+X12+X13
eta2 =~ X21+X22+X23
eta3 =~ X31+X32+X33
eta4 =~ X41+X42+X43
eta5 =~ X51+X52+X53
eta6 =~ X61+X62+X63

# Regressions
eta5 ~ eta1 + eta2 + eta6
eta6 ~ eta3 + eta4 + eta5

# residual covariances
eta5 ~~ eta6
'






