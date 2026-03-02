library(cSEM)
data(ITFlex)

A_itflex <- list(ITComp = ITFlex[, c("ITCOMP1", "ITCOMP2", "ITCOMP3", "ITCOMP4")],
         Modul  = ITFlex[, c("MOD1", "MOD2", "MOD3", "MOD4")],
         ITConn = ITFlex[, c("ITCONN1", "ITCONN2", "ITCONN3", "ITCONN4")],
         ITPers = ITFlex[, c("ITPSF1", "ITPSF2", "ITPSF3", "ITPSF4")])

C_itflex <- matrix(c(0, 0, 0, 0,
             1, 0, 1, 0,
             1, 0, 0, 0,
             1, 1, 1, 0),
           4, 4, byrow = FALSE)

colnames(C_itflex) <- rownames(C_itflex) <- names(A_itflex)


mode_itflex <- rep("formative", 4)
