library(devtools)
load_all()

rm(list = ls())
data(ECSI) ; ECSI = ECSI/10
A = list(CUSTOMER_E = ECSI[, c("CUEX1", "CUEX2", "CUEX3")],
         PERC_QUAL  = ECSI[, c("PERQ1", "PERQ2", "PERQ3", "PERQ4", "PERQ5", "PERQ6", "PERQ7")],
         PERC_VALUE = ECSI[, c("PERV1", "PERV2")],
         CUSTOMER_S = ECSI[, c("CUSA1", "CUSA2", "CUSA3")],
         CUSTOMER_L = ECSI[, c("CUSL1", "CUSL2", "CUSL3")]
)

############
#  svdSEM  #
############

C = matrix(c(0, 0, 0, 0, 0,
             1, 0, 0, 0, 0,
             1, 1, 0, 0, 0,
             1, 1, 1, 0, 0,
             0, 0, 0, 1, 0),
           5, 5, byrow = FALSE)

colnames(C) = rownames(C) = names(A)

svd_fc <- SemFC$new(data=A,
                    relation_matrix = C,
                    mode=rep("reflective", 5),
                    scale = T,
                    estimator = "svd")

svd_fc$fit(infer = T, B = 100)
svd_fc$summary()





C = matrix(c(0, 0, 0, 0, 0,
             1, 0, 0, 0, 0,
             1, 1, 0, 0, 0,
             1, 1, 1, 0, 0,
             0, 0, 0, 1, 0),
           5, 5, byrow = FALSE)

colnames(C) = rownames(C) = names(A)

mode = rep("reflective", 5)
model_fc <- SemFC$new(data=A, relation_matrix = C, mode=mode, scale = T)
model_fc$fit(infer = F)
model_fc$summary()


# Full ECSI
data("ECSI")
ECSI = ECSI/10

L = list(IMAG = ECSI[, 1:5],
         CUEX = ECSI[, 6:8],
         PERQ = ECSI[, 9:15],
         PERV = ECSI[, 16:17],
         CUSA = ECSI[, 18:20],
         CUSCO = ECSI[, 21, drop = F],
         CUSL = ECSI[, 22:24])


C <- matrix(c(0, 1, 0, 0, 1, 0, 1,
              0, 0, 1, 1, 1, 0, 0,
              0, 0, 0, 1, 1, 0, 0,
              0, 0, 0, 0, 1, 0, 0,
              0, 0, 0, 0, 0, 1, 1,
              0, 0, 0, 0, 0, 0, 1,
              0, 0, 0, 0, 0, 0, 0), 7, 7, byrow = TRUE)

colnames(C) = rownames(C) = names(L)

mode = rep("reflective", 7) ; mode[6] = "formative"

model <- SemFC$new(data=L, relation_matrix = C, mode=mode, estimator = "svd")
model$fit(infer = F)
model$fit(infer = T, B = 100)
model$summary()



modelml <- SemFC$new(data=L, relation_matrix = C, mode=mode, estimator = "ml")
modelml$fit(infer = F)
modelml$fit(infer = T)
modelml$summary()