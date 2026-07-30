data(ECSI)
ECSI <- ECSI/10


A <- list(IMAG = ECSI[, 1:5],
          CUEX = ECSI[, 6:8],
          PERQ = ECSI[, 9:15],
          PERV = ECSI[, 16:17],
          CUSA = ECSI[, 18:20],
          COMP = ECSI[, 21, drop = F],
          CUSL = ECSI[, 22:24])


C_ecsi_full <- matrix(c(0, 1, 0, 0, 1, 0, 1,
                        0, 0, 1, 1, 1, 0, 0,
                        0, 0, 0, 1, 1, 0, 0,
                        0, 0, 0, 0, 1, 0, 0,
                        0, 0, 0, 0, 0, 1, 1,
                        0, 0, 0, 0, 0, 0, 1,
                        0, 0, 0, 0, 0, 0, 0), 7, 7, byrow = TRUE)



colnames(C_ecsi_full) <- rownames(C_ecsi_full) <- names(A)

mode_ecsi_full <- rep("reflective", 7) ; mode_ecsi_full[6] <- "formative"

A2 <- data.frame(Reduce("cbind", A))

sem.model.ecsi_full <- '
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