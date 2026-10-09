# Configuration dataset BergamiBagozzi2000
# Source: inst/model/model_bergami.R

create_bergami_config <- function(dataset) {
  A <- list(
    OrgPres = dataset[, c("cei1", "cei2", "cei3", "cei4", "cei5",
                          "cei6", "cei7", "cei8")],
    OrgIden = dataset[, c("ma1", "ma2", "ma3", "ma4", "ma5", "ma6")],
    AffLove = dataset[, c("orgcmt1", "orgcmt2", "orgcmt3", "orgcmt7")],
    AffJoy = dataset[, c("orgcmt5", "orgcmt8")],
    Gender = dataset[, "gender", drop = FALSE]
  )

  C <- matrix(c(0, 0, 0, 0, 0,
                1, 0, 0, 0, 0,
                1, 1, 0, 0, 1,
                1, 1, 0, 0, 1,
                0, 0, 0, 0, 0),
              5, 5, byrow = FALSE)

  colnames(C) <- rownames(C) <- names(A)

  # Mode: reflective pour les 4 premiers, formative pour Gender (selon R/data.R ligne 276)
  mode <- c(rep("reflective", 4), "formative")
  names(mode) <- names(A)


  sem <-'
  # Measurement models
  OrgPres =~ cei1 + cei2 + cei3 + cei4 + cei5 + cei6 + cei7 + cei8
  OrgIden =~ ma1 + ma2 + ma3 + ma4 + ma5 + ma6
  AffLove =~ orgcmt1 + orgcmt2 + orgcmt3 + orgcmt7
  AffJoy  =~ orgcmt5 + orgcmt8
  Gender  <~ gender

  # Structural model
  OrgIden ~ OrgPres
  AffLove ~ OrgPres + OrgIden + Gender
  AffJoy  ~ OrgPres + OrgIden + Gender


  #covariances
  AffLove ~~ 0*AffJoy


  '

  return(list(data = A, relation_matrix = C, mode = mode, sem = sem))
}
