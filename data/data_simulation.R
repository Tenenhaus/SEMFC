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

R = rbind(cbind(PHI, PHI%*%t(GAMMA)%*%t(solve(diag(NROW(BETA))-BETA))),
          cbind(solve(diag(NROW(BETA))-BETA)%*%GAMMA%*%PHI,
                solve(diag(NROW(BETA))-BETA)%*%(GAMMA%*%PHI%*%t(GAMMA)+PSI)
                %*%t(solve(diag(NROW(BETA))-BETA)))
)


R2_1 = 1-PSI[1, 1] ; R2_2 = 1-PSI[2, 2] ; R2 = c(R2_1, R2_2)

SIGMA11 = SIGMA22 = SIGMA33 = SIGMA44 = matrix(c(1, .3, .4,
                                                 .3, 1, .5,
                                                 .4, .5, 1), 3, 3)

SIGMA55 = SIGMA66 = matrix(c(1, .49, .49,
                             .49, 1, .49,
                             .49, .49, 1), 3, 3)

w_exo_1 = w_exo_2 = rep(1, 3)/drop(sqrt(t(rep(1, 3))%*%SIGMA11%*%rep(1, 3)))
w_exo_3 = w_exo_4 = (1:3)/drop(sqrt(t(1:3)%*%SIGMA11%*%(1:3)))

l1 = l2 = SIGMA11%*%w_exo_1
l3 = l4 = SIGMA11%*%w_exo_3
l5 = l6 = rep(.7, 3)
LAMBDA = Matrix::bdiag(list(l1, l2, l3, l4, l5, l6))
lambda = list(l1, l2, l3, l4, l5, l6)
omega = list(w_exo_1, w_exo_2, w_exo_3, w_exo_4)

SIGMA = Matrix::bdiag(list(SIGMA11, SIGMA22, SIGMA33,
                   SIGMA44, SIGMA55, SIGMA66))

SIGMA[1:3, 4:6]   = R[1,2]*l1%*%t(l2);SIGMA[4:6, 1:3] = t(SIGMA[1:3, 4:6])
SIGMA[1:3, 7:9]   = R[1,3]*l1%*%t(l3);SIGMA[7:9, 1:3] = t(SIGMA[1:3, 7:9])
SIGMA[1:3, 10:12] = R[1,4]*l1%*%t(l4);SIGMA[10:12, 1:3] = t(SIGMA[1:3, 10:12])
SIGMA[1:3, 13:15] = R[1,5]*l1%*%t(l5);SIGMA[13:15, 1:3] = t(SIGMA[1:3, 13:15])
SIGMA[1:3, 16:18] = R[1,6]*l1%*%t(l6);SIGMA[16:18, 1:3] = t(SIGMA[1:3, 16:18])

SIGMA[4:6, 7:9] = R[2,3]*l2%*%t(l3);SIGMA[7:9, 4:6] = t(SIGMA[4:6, 7:9])
SIGMA[4:6, 10:12] = R[2,4]*l2%*%t(l4);SIGMA[10:12, 4:6] = t(SIGMA[4:6, 10:12])
SIGMA[4:6, 13:15] = R[2,5]*l2%*%t(l5);SIGMA[13:15, 4:6] = t(SIGMA[4:6, 13:15])
SIGMA[4:6, 16:18] = R[2,6]*l2%*%t(l6);SIGMA[16:18, 4:6] = t(SIGMA[4:6, 16:18])

SIGMA[7:9, 10:12] = R[3,4]*l3%*%t(l4);SIGMA[10:12, 7:9] = t(SIGMA[7:9, 10:12])
SIGMA[7:9, 13:15] = R[3,5]*l3%*%t(l5);SIGMA[13:15, 7:9] = t(SIGMA[7:9, 13:15])
SIGMA[7:9, 16:18] = R[3,6]*l3%*%t(l6);SIGMA[16:18, 7:9] = t(SIGMA[7:9, 16:18])

SIGMA[10:12, 13:15] = R[4,5]*l4%*%t(l5);SIGMA[13:15, 10:12] = t(SIGMA[10:12, 13:15])
SIGMA[10:12, 16:18] = R[4,6]*l4%*%t(l6);SIGMA[16:18, 10:12] = t(SIGMA[10:12, 16:18])

SIGMA[13:15, 16:18] = R[5,6]*l5%*%t(l6);SIGMA[16:18, 13:15] = t(SIGMA[13:15, 16:18])

true_param_with_S = c(l1, l2 , l3, l4, l5, l6,
                      R[1:4, 1:4][upper.tri(R[1:4, 1:4])],
                      GAMMA[1, 1:2], GAMMA[2, 3:4],
                      BETA[1, 2], BETA[2, 1],
                      R[5, 6],
                      SIGMA11[upper.tri(SIGMA11, diag = TRUE)],
                      SIGMA22[upper.tri(SIGMA22, diag = TRUE)],
                      SIGMA33[upper.tri(SIGMA33, diag = TRUE)],
                      SIGMA44[upper.tri(SIGMA44, diag = TRUE)],
                      1 - l5^2,
                      1 - l6^2
)
