


parameters_svd <- function(lambda,
                           P_EXO,
                           G,
                           B,
                           P_ENDO,
                           residual_variance,
                           S_composites,
                           mode){

  J <- length(mode)

  S_composites_upper <- lapply(S_composites, function(matrix) matrix[upper.tri(matrix, diag = TRUE)])

  # empirical covariance for composite blocks or residual_variance for reflective
  diag_jj <- vector("list", J)
  diag_jj[mode == "formative"] <- S_composites_upper
  diag_jj[mode != "formative"] <- residual_variance


  theta_vect <-
    Reduce("c",
      c(Reduce("c", lambda),
        P_EXO[upper.tri(P_EXO)],
        apply(G, 1, function(row) row[row != 0]),
        apply(B, 1, function(row) row[row != 0]),
        P_ENDO[upper.tri(P_ENDO)],
        Reduce("c", diag_jj)

      )
    )

  return(unname(theta_vect))
}