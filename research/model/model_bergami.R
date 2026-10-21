library(cSEM)
data("BergamiBagozzi2000")


A_bergami <- list(OrgPres = BergamiBagozzi2000[, c("cei1", "cei2", "cei3", "cei4", "cei5",
                          "cei6", "cei7", "cei8")],
         OrgIden = BergamiBagozzi2000[, c("ma1", "ma2", "ma3", "ma4", "ma5", "ma6")],
         AffLove = BergamiBagozzi2000[, c("orgcmt1", "orgcmt2", "orgcmt3", "orgcmt7")],
         AffJoy  = BergamiBagozzi2000[, c("orgcmt5", "orgcmt8")],
         Gender  = BergamiBagozzi2000[, "gender", drop = FALSE]
         )

C_bergami <- matrix(c(0, 0, 0, 0, 0,
             1, 0, 0, 0, 0,
             1, 1, 0, 0, 1,
             1, 1, 0, 0, 1,
             0, 0, 0, 0, 0),
           5, 5, byrow = FALSE)

colnames(C_bergami) <- rownames(C_bergami) <- names(A_bergami)


mode_bergami <- c(rep("reflective", 4), "formative")



sem.model.bergami <-'
# Measurement models
OrgPres =~ cei1 + cei2 + cei3 + cei4 + cei5 + cei6 + cei7 + cei8
OrgIden =~ ma1 + ma2 + ma3 + ma4 + ma5 + ma6
AffLove =~ orgcmt1 + orgcmt2 + orgcmt3 + orgcmt7
AffJoy  =~ orgcmt5 + orgcmt8
Gender  =~ gender

# Structural model
OrgIden ~ OrgPres
AffLove ~ OrgPres + OrgIden + Gender
AffJoy  ~ OrgPres + OrgIden + Gender


#covariances
AffLove ~~ 0*AffJoy


'




