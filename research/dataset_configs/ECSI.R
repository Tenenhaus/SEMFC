# Configuration dataset ECSI
# Retourne une liste avec:
#   - data: liste des matrices de variables manifestes
#   - relation_matrix: matrice de relations entre construits

create_ecsi_config <- function(dataset) {
  # Normaliser ECSI
  dataset <- dataset / 10

  A <- list(
    CUSTOMER_E = dataset[, c("CUEX1", "CUEX2", "CUEX3")],
    PERC_QUAL = dataset[, c("PERQ1", "PERQ2", "PERQ3", "PERQ4", "PERQ5", "PERQ6", "PERQ7")],
    PERC_VALUE = dataset[, c("PERV1", "PERV2")],
    CUSTOMER_S = dataset[, c("CUSA1", "CUSA2", "CUSA3")],
    CUSTOMER_L = dataset[, c("CUSL1", "CUSL2", "CUSL3")]
  )

  C <- matrix(c(0, 0, 0, 0, 0,
                 1, 0, 0, 0, 0,
                 1, 1, 0, 0, 0,
                 1, 1, 1, 0, 0,
                 0, 0, 0, 1, 0), 5, 5, byrow = FALSE)
  colnames(C) <- rownames(C) <- names(A)

  # Mode de chaque bloc latent (reflective)
  mode <- rep("reflective", length(A))
  names(mode) <- names(A)


  sem <-  '
  # latent variable definitions
      CUSTOMER_E =~ CUEX1+CUEX2+CUEX3
      PERC_QUAL =~ PERQ1+PERQ2+PERQ3+PERQ4+PERQ5+PERQ6+PERQ7
      PERC_VALUE =~ PERV1+PERV2
      CUSTOMER_S =~ CUSA1+CUSA2+CUSA3
      CUSTOMER_L =~ CUSL1+CUSL2+CUSL3

      # Regressions
      PERC_QUAL ~ CUSTOMER_E
      PERC_VALUE ~ CUSTOMER_E + PERC_QUAL
      CUSTOMER_S ~ CUSTOMER_E + PERC_QUAL + PERC_VALUE
      CUSTOMER_L ~ CUSTOMER_S'

  return(list(data = A, relation_matrix = C, mode = mode, sem = sem))
}


