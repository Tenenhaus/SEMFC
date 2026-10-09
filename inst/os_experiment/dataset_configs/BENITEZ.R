# Configuration dataset Benitezetal2020
# Source: inst/model/model_benitez.R

create_benitez_config <- function(dataset) {
  A = list(SEXB = dataset[, c("SEXB1", "SEXB2", "SEXB3", "SEXB4")],
           SEMB  = dataset[, c("SEMB1", "SEMB2", "SEMB3", "SEMB4")],
           SMC = dataset[, c("SMC1","SMC2", "SMC3", "SMC4")],
           BPP = dataset[, c("BPP1", "BPP2", "BPP3", "BPP4", "BPP5")],
           FS = dataset[, "FirmSize", drop = F],
           Ind = dataset[, c("Industry1", "Industry2", "Industry3")])


  C = matrix(c(0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0,
               1, 1, 0, 0, 0, 0,
               0, 0, 1, 0, 1, 1,
               0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0),
             6, 6, byrow = FALSE)


  colnames(C) <- rownames(C) <- names(A)

  # Mode: reflective pour les 4 premiers, formative pour Gender (selon R/data.R ligne 276)
  mode = c(rep("reflective", 2), rep("formative", 4))
  names(mode) <- names(A)


  sem <-'
  # Reflective measurement models
  SEXB =~ SEXB1 + SEXB2 + SEXB3 +SEXB4
  SEMB =~ SEMB1 + SEMB2 + SEMB3 + SEMB4

  # Composite models
  SMC <~ SMC1 + SMC2 + SMC3 + SMC4
  BPP <~ BPP1 + BPP2 + BPP3 + BPP4 + BPP5

  # Control variables
  FS <~ FirmSize
  Ind <~ Industry1 + Industry2 + Industry3

  # Structural model
  SMC ~ SEXB + SEMB
  BPP ~ SMC + Ind + FS
  '

  return(list(data = A, relation_matrix = C, mode = mode, sem = sem))
}
