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

mode_ecsi <- rep("reflective", length(A))


A2 <- data.frame(Reduce("cbind", A))

sem.model.ecsi <-  '
# latent variable definitions
    CUSTOMER_E =~ CUEX1+CUEX2+CUEX3
    PERC_QUAL =~ PERQ1+PERQ2+PERQ3+PERQ4+PERQ5+PERQ6+PERQ7
    PERC_VALUE =~ PERV1+PERV2
    CUSTOMER_S =~ CUSA1+CUSA2+CUSA3
    CUSTOMER_L =~ CUSL1+CUSL2+CUSL3

    # Regressions
    PERC_QUAL ~ CUSTOMER_E
    PERC_VALUE ~ CUSTOMER_E + PERC_QUAL
    CUSTOMER_S ~ CUSTOMER_E + PERC_QUAL + PERC_VALUE
    CUSTOMER_L ~ CUSTOMER_S'