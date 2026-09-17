# Configuration dataset LancelotMiltgenetal2016
# Source: inst/examples/lancelot.R

create_lancelot_config <- function(dataset) {
  A <- list(
    Trust = dataset[, grep("trust", colnames(dataset), ignore.case = TRUE)],
    PrCon = dataset[, grep("priv", colnames(dataset), ignore.case = TRUE)],
    Risk = dataset[, grep("risk", colnames(dataset), ignore.case = TRUE)],
    Intent = dataset[, grep("intent", colnames(dataset), ignore.case = TRUE)]
  )

  C <- matrix(c(0, 1, 0, 0,
                0, 0, 0, 0,
                1, 1, 0, 0,
                1, 1, 1, 0),
              4, 4, byrow = FALSE)

  colnames(C) <- rownames(C) <- names(A)

  # Mode: reflective (selon R/data.R ligne 483)
  mode <- rep("reflective", length(A))
  names(mode) <- names(A)

  return(list(data = A, relation_matrix = C, mode = mode))
}


