# Configuration dataset ITFlex
# Source: inst/model/model_itflex.R

create_itflex_config <- function(dataset) {
  A <- list(
    ITComp = dataset[, c("ITCOMP1", "ITCOMP2", "ITCOMP3", "ITCOMP4")],
    Modul = dataset[, c("MOD1", "MOD2", "MOD3", "MOD4")],
    ITConn = dataset[, c("ITCONN1", "ITCONN2", "ITCONN3", "ITCONN4")],
    ITPers = dataset[, c("ITPSF1", "ITPSF2", "ITPSF3", "ITPSF4")]
  )

  C <- matrix(c(0, 0, 0, 0,
                1, 0, 1, 0,
                1, 0, 0, 0,
                1, 1, 1, 0),
              4, 4, byrow = FALSE)

  colnames(C) <- rownames(C) <- names(A)

  # Mode: formative (selon R/data.R ligne 385)
  mode <- rep("formative", length(A))
  names(mode) <- names(A)

  sem <- '
  ITComp <~ ITCOMP1 + ITCOMP2 + ITCOMP3 + ITCOMP4
  Modul  <~ MOD1 + MOD2 + MOD3 + MOD4
  ITConn <~ ITCONN1 + ITCONN2 + ITCONN3 + ITCONN4
  ITPers <~ ITPSF1 + ITPSF2 + ITPSF3 + ITPSF4

  Modul  ~ ITComp + ITConn
  ITConn ~ ITComp
  ITPers ~ ITComp + Modul + ITConn
'


  return(list(data = A, relation_matrix = C, mode = mode, sem = sem))
}
