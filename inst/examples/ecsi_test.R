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


library(lavaan)

lavaan_ecsi <- '
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

lavaan_ml <- sem(lavaan_ecsi,
                  data=ECSI,
                  estimator = "ML",
                  likelihood="wishart")

summary(lavaan_ml, standardized=TRUE, fit.measures=TRUE)
estimate = parameterEstimates(lavaan_ml, standardized = TRUE)



std_all_ml = unlist(modelml$estimate$std_lambda)
std_all_lavaan = estimate[1:24,11]



lambda_comparaison_ecsi = cbind(estimate[1:24,1:3], round(std_all_ml,3), round(std_all_lavaan,3))
print('lambda')
print(lambda_comparaison_ecsi)



g =  unlist(modelml$estimate$gamma)
b = unlist(modelml$estimate$beta)
bg_ml = c(g[1], b[2,1], b[3,1],b[3,2], g[4],b[4,1],b[4,2],b[4,3], b[5,4], g[6], b[6,4], b[6,5])

bg_lavaan = estimate[25:36,11]

beta_gama_comparaison_ecsi = cbind(estimate[25:36,1:3], round(bg_ml,3), round(bg_lavaan,3))
print('beta et gamma')
print(beta_gama_comparaison_ecsi)