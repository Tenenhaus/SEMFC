# correction of the RGCCA estimates

# correction = function(A, a, bias = FALSE){
#   S = cov2(A, bias = bias)
#   diag(S) = 0
#   d = drop(sqrt((t(a)%*%S%*%a)/(1-sum(a^4))))
#   return(list(d = d, n = 0))
# }

correction = function(A, a, mode = "reflective", bias = FALSE){
  
  if(mode == "reflective"){
    S = cov2(A, bias = bias)
    W = a%*%t(a);diag(W) = 0
    WS = W*S
    d = ifelse(any(WS>0), sqrt(sum(WS[WS>0])/sum(W[WS>0]^2)), 1)
    #d =  sqrt(sum(WS)/sum(W^2))
    n = length(which(WS<0))
  }
  
  if(mode == "formative"){
    #S = cov2(A, bias = bias)
    #eig = eigen(S)
    #Sinv_sqrt  = eig$vectors %*% diag(eig$values^(-1/2)) %*% t(eig$vectors)
    #d = 1/norm(Sinv_sqrt%*%a, type = "2")
    #n = sum(eig$values == 0)
    d = 1/drop(sqrt(t(a)%*%solve(cov2(A, bias = bias))%*%a))
    n = 0
  }
  
  return(list(d = d, n = n))
}
  
  
