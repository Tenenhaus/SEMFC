source('data/data_simulation_mixed.R')


set.seed(20091979)
N <- 300
library(mvtnorm)
X <- mvrnorm(N, rep(0, ncol(SIGMA)), SIGMA, empirical = TRUE)
colnames(X) <- paste0("X",
                      rep(seq_along(lengths(lambda)), times = lengths(lambda)),
                      unlist(lapply(lengths(lambda), seq_len)))


Y <- list(LV1 = X[, range_index[[1]]], LV2 = X[, range_index[[2]]], LV3 = X[, range_index[[3]]],
         LV4 = X[,range_index[[4]]], LV5 = X[, range_index[[5]]], LV6 = X[, range_index[[6]]])





C <- matrix(c(0, 0, 0, 0, 1, 0,
             0, 0, 0, 0, 1, 0,
             0, 0, 0, 0, 0, 1,
             0, 0, 0, 0, 0, 1,
             0, 0, 0, 0, 0, 1,
             0, 0, 0, 0, 1, 0), 6, 6, byrow = TRUE)

colnames(C) <- rownames(C) <- names(Y)

mode <- c(rep("formative", 4), rep("reflective", 2))


X_2 <- mvrnorm(N,  rep(0, ncol(SIGMA)), SIGMA, empirical = FALSE)
colnames(X_2) <- colnames(X)


Y_2 <- list(LV1 = X_2[, range_index[[1]]], LV2 = X_2[, range_index[[2]]], LV3 = X_2[, range_index[[3]]],
         LV4 = X_2[,range_index[[4]]], LV5 = X_2[, range_index[[5]]], LV6 = X_2[, range_index[[6]]])




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


