# restrictions for the minimization of the Loglikelihood function
heq1 <- function(x, S) {
  
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

heq0 <- function(x, S) {
  
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
  
  R = rbind(cbind(P_EXO, P_EXO%*%t(G)%*%t(solve(diag(2) - B))),
            cbind(solve(diag(2) - B)%*%G%*%P_EXO, P_ENDO))
  
  h <- c(rep(0,5))
  h[1] <- t(l1)%*%solve(S1)%*%l1 - 1
  h[2] <- t(l2)%*%solve(S2)%*%l2 - 1
  h[3] <- t(l3)%*%solve(S3)%*%l3 - 1
  h[4] <- t(l4)%*%solve(S4)%*%l4 - 1
  h[5] <- x[27] - x[26]
  
  return(h)
}



