

# pour bloc formatif
reliability <- function (metric, lambdas, residual_variances){
  if (metric=='Dillon'){
    # Calculer la somme des loadings et des résidus
    sum_loadings <- sapply(lambdas, function(x) sum(as.numeric(x)))
    sum_residuals <- sapply(residual_variances, function(x) sum(as.numeric(x)))

    # Calculer numérateur et dénominateur
    numerator <- sum_loadings^2
    denominator <- numerator + sum_residuals

    # Calculer rho pour chaque bloc
    results <- numerator / denominator


    names(results) <- names(lambdas)
    return(results)

  }


}

