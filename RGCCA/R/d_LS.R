#d_LS function
# Euclidean Distance

d_LS = function(A, B){
  A = as.matrix(A)
  B = as.matrix(B)
  d_LS = sqrt(0.5*sum(diag((B-A)%*%(B-A))))
  return(d_LS = d_LS)
}