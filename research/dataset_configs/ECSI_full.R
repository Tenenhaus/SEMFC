# Configuration dataset ECSI_full
# Même données que ECSI mais avec une structure de modèle étendue

create_ecsi_full_config <- function(dataset) {
  # Normaliser ECSI
  dataset <- dataset / 10

  A <- list(
    IMAG = dataset[, c("IMAG1", "IMAG2", "IMAG3", "IMAG4", "IMAG5")],
    CUEX = dataset[, c("CUEX1", "CUEX2", "CUEX3")],
    PERQ = dataset[, c("PERQ1", "PERQ2", "PERQ3", "PERQ4", "PERQ5", "PERQ6", "PERQ7")],
    PERV = dataset[, c("PERV1", "PERV2")],
    CUSA = dataset[, c("CUSA1", "CUSA2", "CUSA3")],
    COMP = dataset[, c("CUSCO"), drop = FALSE],
    CUSL = dataset[, c("CUSL1", "CUSL2", "CUSL3")]
  )

  C <- matrix(c(0, 1, 0, 0, 1, 0, 1,
                0, 0, 1, 1, 1, 0, 0,
                0, 0, 0, 1, 1, 0, 0,
                0, 0, 0, 0, 1, 0, 0,
                0, 0, 0, 0, 0, 1, 1,
                0, 0, 0, 0, 0, 0, 1,
                0, 0, 0, 0, 0, 0, 0), 7, 7, byrow = TRUE)

  colnames(C) <- rownames(C) <- names(A)

  # Mode: mostly reflective except CUSCO which is formative
  mode <- rep("reflective", 7)
  mode[6] <- "formative"
  names(mode) <- names(A)

  sem <- '
    # measurement model
    IMAG  =~ IMAG1 + IMAG2 + IMAG3 + IMAG4 + IMAG5
    CUEX  =~ CUEX1 + CUEX2 + CUEX3
    PERQ  =~ PERQ1 + PERQ2 + PERQ3 + PERQ4 + PERQ5 + PERQ6 + PERQ7
    PERV  =~ PERV1 + PERV2
    CUSA  =~ CUSA1 + CUSA2 + CUSA3
    COMP =~ CUSCO
    CUSL  =~ CUSL1 + CUSL2 + CUSL3

    # structural model
    CUEX ~ IMAG
    PERQ ~ CUEX
    PERV ~ CUEX + PERQ
    CUSA ~ IMAG + CUEX + PERQ + PERV
    COMP ~ CUSA
    CUSL ~ IMAG + CUSA + COMP
    '

  return(list(data = A, relation_matrix = C, mode = mode, sem = sem))
}

