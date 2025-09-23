######################################
####       data simulation        ####
#### true value of the parameters ####
######################################
library(Matrix)

BETA = matrix(c(0, 0.25,
                0.5,  0), 2, 2, byrow = TRUE)

GAMMA = matrix(c(-0.30, 0.5, 0, 0,
                 0,  0,  0.5, 0.25), 2, 4, byrow = TRUE)

R22 = matrix(c(1, sqrt(1/2), sqrt(1/2), 1), 2, 2)

PHI = matrix(c( 1, .5,  .5, .5,
                .5,  1,  .5, .5,
                .5, .5,   1, .5,
                .5, .5,  .5,  1), 4, 4, byrow = TRUE)

PSI = (diag(2)-BETA)%*%R22%*%t(diag(2)-BETA) - GAMMA%*%PHI%*%t(GAMMA)

R = rbind(
  cbind(PHI,
        PHI%*%t(GAMMA)%*%t(solve(diag(NROW(BETA))-BETA))),
  cbind(solve(diag(NROW(BETA))-BETA)%*%GAMMA%*%PHI,
        solve(diag(NROW(BETA))-BETA)%*%(GAMMA%*%PHI%*%t(GAMMA)+PSI) %*%t(solve(diag(NROW(BETA))-BETA)))
)


R2 = 1- diag(PSI)



l1 = c(rep(0.9, 50) ,rep(0, 100), rep(0.2,30))
#l1 = c(rep(0.9, 10) ,rep(0, 10))
# l1 = rep(.7, 3)
# l1 = c(rep(0.9, 3), 0.2, rep(0, 5))
l2 = rep(.7, 3)
l3 = l4 = rep(.7, 3)
l5 = l6 = rep(.7, 3)

lambda = list(l1, l2, l3, l4, l5, l6)
LAMBDA = Matrix::bdiag(lambda)

J = length(lambda)

SIGMA_diag = list()
for (j in 1:J){
  SIGMAjj = lambda[[j]]%*% t(lambda[[j]])
  diag(SIGMAjj) = rep(1, ncol(SIGMAjj))
  SIGMA_diag <- append(SIGMA_diag, list(SIGMAjj))

}



# w_exo_1 = w_exo_2 = rep(1, 3)/drop(sqrt(t(rep(1, 3))%*%SIGMA11%*%rep(1, 3)))
# w_exo_3 = w_exo_4 = (1:3)/drop(sqrt(t(1:3)%*%SIGMA11%*%(1:3)))

# omega = list(w_exo_1, w_exo_2, w_exo_3, w_exo_4)

SIGMA = Matrix::bdiag(SIGMA_diag)

index_end = cumsum(lengths(lambda))
index_start = c(1, index_end[-length(index_end)] + 1)
range_index <- lapply(1:6, function(i) index_start[i]:index_end[i])

for (j in 1:(J-1)){

  for (i in (j+1):J){
    range_row = range_index[[j]]
    range_col = range_index[[i]]

    lj = lambda[[j]]
    li = lambda[[i]]

    SIGMA[range_row, range_col] = R[j,i]*lj%*%t(li)
    SIGMA[range_col, range_row] = t(SIGMA[range_row, range_col] )


  }


}



true_param_with_S = c(l1, l2 , l3, l4, l5, l6,
                      R[1:4, 1:4][upper.tri(R[1:4, 1:4])],
                      GAMMA[1, 1:2], GAMMA[2, 3:4],
                      BETA[1, 2], BETA[2, 1],
                      R[5, 6],
                      1 - l1^2,
                      1 - l2^2,
                      1 - l3^2,
                      1 - l4^2,
                      1 - l5^2,
                      1 - l6^2
)
