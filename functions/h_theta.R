###################################
#      return h(\theta) for       #
#       hessian computation       #
###################################

h_theta <- function(x, block = 1){
  #loadings
  L = list(l1 = x[1:3], l2 = x[4:6], 
           l3 = x[7:9], l4 = x[10:12])
  
  S = list( 
    #cov between MVs
    S1 = matrix(c(x[32], x[33], x[35],
                  x[33], x[34], x[36],
                  x[35], x[36], x[37]), 3, 3),
    #cov between MVs
    S2 = matrix(c(x[38], x[39], x[41],
                  x[39], x[40], x[42],
                  x[41], x[42], x[43]), 3, 3),
    #cov between MVs
    S3 = matrix(c(x[44], x[45], x[47],
                  x[45], x[46], x[48],
                  x[47], x[48], x[49]), 3, 3),
    #cov between MVs
    S4 = matrix(c(x[50], x[51], x[53],
                  x[51], x[52], x[54],
                  x[53], x[54], x[55]), 3, 3)
  )
  
  h = drop(L[[block]]%*%solve(S[[block]])%*%L[[block]] - 1)
  
}