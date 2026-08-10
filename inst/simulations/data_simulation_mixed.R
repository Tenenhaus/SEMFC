######################################
####       data simulation        ####
#### true value of the parameters ####
######################################
library(Matrix)

mode <- c(rep("formative", 4), rep("reflective", 2))
J <- length(mode)

BETA <- matrix(c(0, 0.25,
                0.5,  0), 2, 2, byrow = TRUE)

GAMMA <- matrix(c(-0.30, 0.5, 0, 0,
                 0,  0,  0.5, 0.25), 2, 4, byrow = TRUE)

R22 <- matrix(c(1, sqrt(1/2),
               sqrt(1/2), 1), 2, 2)

PHI <- matrix(c( 1, .5,  .5, .5,
                .5,  1,  .5, .5,
                .5, .5,   1, .5,
                .5, .5,  .5,  1), 4, 4, byrow = TRUE)

PSI <- (diag(2)-BETA)%*%R22%*%t(diag(2)-BETA) - GAMMA%*%PHI%*%t(GAMMA)

R <- rbind(cbind(PHI, PHI%*%t(GAMMA)%*%t(solve(diag(NROW(BETA))-BETA))),
          cbind(solve(diag(NROW(BETA))-BETA)%*%GAMMA%*%PHI,
                solve(diag(NROW(BETA))-BETA)%*%(GAMMA%*%PHI%*%t(GAMMA)+PSI)
                %*%t(solve(diag(NROW(BETA))-BETA)))
)


# R2_1 <- 1-PSI[1, 1] ; R2_2 <- 1-PSI[2, 2]
R2 <- 1- diag(PSI)


# q <- 3
rho <- 0.4

SIGMA_form <- matrix(rho, q, q)
diag(SIGMA_form) <- 1

SIGMA_form[1:3, 1:3] <- matrix(
  c(
    1,   0.3, 0.4,
    0.3, 1,   0.5,
    0.4, 0.5, 1
  ),
  nrow = 3,
  byrow = TRUE
)

SIGMA11 <- SIGMA_form
SIGMA22 <- SIGMA_form
SIGMA33 <- SIGMA_form
SIGMA44 <- SIGMA_form



SIGMA_ref <- matrix(0.49, q, q)
diag(SIGMA_ref) <- 1

SIGMA55 <- SIGMA_ref
SIGMA66 <- SIGMA_ref






w_exo_1 <- rep(1, q)
w_exo_1  <- w_exo_1/drop(sqrt(t(w_exo_1)%*%SIGMA11%*%w_exo_1))

w_exo_2 <- rep(1, q)
w_exo_2  <- w_exo_2/drop(sqrt(t(w_exo_2)%*%SIGMA22%*%w_exo_2))


w_exo_3 <- seq(1, 3, length.out = q)
w_exo_3 <- w_exo_4 <- w_exo_3/drop(sqrt(t(w_exo_3)%*%SIGMA22%*%w_exo_3))

l1 <- SIGMA11%*%w_exo_1
l2 <- SIGMA22%*%w_exo_2
l3 <- l4 <- SIGMA22%*%w_exo_3
l5 <- l6 <- rep(.7, q)


lambda <- list(l1, l2, l3, l4, l5, l6)
LAMBDA <- Matrix::bdiag(lambda)


omega <- list(w_exo_1, w_exo_2, w_exo_3, w_exo_4)

SIGMA <- Matrix::bdiag(list(SIGMA11, SIGMA22, SIGMA33,
                   SIGMA44, SIGMA55, SIGMA66))

index_end <- cumsum(lengths(lambda))
index_start <- c(1, index_end[-length(index_end)] + 1)
range_index <- lapply(1:6, function(i) index_start[i]:index_end[i])

for (j in 1:(J-1)){

  for (i in (j+1):J){
    range_row <- range_index[[j]]
    range_col <- range_index[[i]]

    lj <- lambda[[j]]
    li <- lambda[[i]]

    SIGMA[range_row, range_col] <- R[j,i]*lj%*%t(li)
    SIGMA[range_col, range_row] <- t(SIGMA[range_row, range_col] )


  }


}

true_param_with_S <- c(l1, l2 , l3, l4, l5, l6,
                      R[1:4, 1:4][upper.tri(R[1:4, 1:4])],
                      GAMMA[1, 1:2], GAMMA[2, 3:4],
                      BETA[2, 1], BETA[1, 2],
                      R[5, 6],
                      SIGMA11[lower.tri(SIGMA11, diag = TRUE)],
                      SIGMA22[lower.tri(SIGMA22, diag = TRUE)],
                      SIGMA33[lower.tri(SIGMA33, diag = TRUE)],
                      SIGMA44[lower.tri(SIGMA44, diag = TRUE)],
                      1 - l5^2,
                      1 - l6^2
)


