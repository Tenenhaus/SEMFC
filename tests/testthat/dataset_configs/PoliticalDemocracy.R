# Configuration dataset PoliticalDemocracy
# Source: R/data.R documentation

create_politicaldemocracy_config <- function(dataset) {
  A <- list(
    ind60 = dataset[, c("x1", "x2", "x3")],
    dem60 = dataset[, c("y1", "y2", "y3", "y4")],
    dem65 = dataset[, c("y5", "y6", "y7", "y8")]
  )

  C <- matrix(c(0, 0, 0,
                1, 0, 0,
                1, 1, 0),
              3, 3, byrow = FALSE)

  colnames(C) <- rownames(C) <- names(A)

  # Mode: reflective (selon R/data.R ligne 565)
  mode <- rep("reflective", length(A))
  names(mode) <- names(A)

  return(list(data = A, relation_matrix = C, mode = mode))
}



