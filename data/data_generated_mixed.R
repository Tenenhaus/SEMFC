source('data/data_simulation_mixed.R')


set.seed(20091979)
N <- 300
library(mvtnorm)
X <- mvrnorm(N, rep(0, 18), SIGMA, empirical = TRUE)
colnames(X) <- paste("X", rep(1:6, each = 3), rep(1:3, 6), sep ="")


Y <- list(X1 = X[, 1:3], X2 = X[, 4:6], X3 = X[, 7:9],
         X4 = X[, 10:12], X5 = X[, 13:15], X6 = X[, 16:18])



C <- matrix(c(0, 0, 0, 0, 1, 0,
             0, 0, 0, 0, 1, 0,
             0, 0, 0, 0, 0, 1,
             0, 0, 0, 0, 0, 1,
             0, 0, 0, 0, 0, 1,
             0, 0, 0, 0, 1, 0), 6, 6, byrow = TRUE)

colnames(C) <- rownames(C) <- names(Y)

mode <- c(rep("formative", 4), rep("reflective", 2))


X_2 <- mvrnorm(N, rep(0, 18), SIGMA, empirical = FALSE)
colnames(X_2) <- paste("X", rep(1:6, each = 3), rep(1:3, 6), sep ="")


Y_2 <- list(LV1 = X_2[, 1:3], LV2 = X_2[, 4:6], LV3 = X_2[, 7:9],
         LV4 = X_2[, 10:12], LV5 = X_2[, 13:15], LV6 = X_2[, 16:18])