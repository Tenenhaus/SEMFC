# Configuration dataset Russett
# Source: R/data.R documentation

create_russett_config <- function(dataset) {
  A <- list(
    AgrIneq = dataset[, c("gini", "farm", "rent")],
    IndDev = dataset[, c("gnpr", "labo")],
    PolInst = dataset[, c("inst", "ecks", "death", "demostab", "dictator")]
  )

  C <- matrix(c(0, 0, 0,
                0, 0, 0,
                1, 1, 0),
              3, 3, byrow = FALSE)

  colnames(C) <- rownames(C) <- names(A)

  # Mode: formative (selon R/data.R ligne 178)
  mode <- rep("formative", length(A))
  names(mode) <- names(A)



  sem <- '
  # 1. Modèle de mesure (Variables latentes)
  AgrIneq <~ gini + farm + rent
  IndDev  <~ gnpr + labo
  PolInst <~ inst + ecks + death + demostab + dictator

  # 2. Modèle structurel (Régressions)
  PolInst ~ AgrIneq + IndDev
'

  return(list(data = A, relation_matrix = C, mode = mode, sem = sem))
}



