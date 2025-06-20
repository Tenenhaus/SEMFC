library(readxl)
ECSI <- as.data.frame(read_excel("data/mobil.xls"))/10
A <- list(CUSTOMER_E = ECSI[, c("CUEX1", "CUEX2", "CUEX3")],
         PERC_QUAL  = ECSI[, c("PERQ1", "PERQ2", "PERQ3", "PERQ4",
                               "PERQ5", "PERQ6", "PERQ7")],
         PERC_VALUE = ECSI[, c("PERV1", "PERV2")],
         CUSTOMER_S = ECSI[, c("CUSA1", "CUSA2", "CUSA3")],
         CUSTOMER_L = ECSI[, c("CUSL1", "CUSL2", "CUSL3")])

C_ecsi <- matrix(c(0, 1, 1, 1, 0,
              0, 0, 1, 1, 0,
              0, 0, 0, 1, 0,
              0, 0, 0, 0, 1,
              0, 0, 0, 0, 0), 5, 5, byrow = TRUE)


colnames(C_ecsi) <- rownames(C_ecsi) <- names(A)

mode_ecsi <- rep("reflective", 5)


A2 <- data.frame(Reduce("cbind", A))

sem.model.ecsi <-  '
# latent variable definitions
    eta1 =~ CUEX1+CUEX2+CUEX3
    eta2 =~ PERQ1+PERQ2+PERQ3+PERQ4+PERQ5+PERQ6+PERQ7
    eta3 =~ PERV1+PERV2
    eta4 =~ CUSA1+CUSA2+CUSA3
    eta5 =~ CUSL1+CUSL2+CUSL3

    # Regressions
    eta2 ~ eta1
    eta3 ~ eta1 + eta2
    eta4 ~ eta1 + eta2 + eta3
    eta5 ~ eta4'