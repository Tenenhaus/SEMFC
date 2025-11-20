
# install.packages("remotes")
remotes::install_github("yrosseel/lavaan")

library(lavaan)

source('R/SEMFC/sem_f_c.R')




source('data/data_generated_mixed.R')



sem.model <-  '
# latent variable definitions
eta5 =~ X51+X52+X53
eta6 =~ X61+X62+X63

# Composite model
eta1 <~ X11 + X12 + X13
eta2 <~ X21 + X22 + X23
eta3 <~ X31 + X32 + X33
eta4 <~ X41 + X42 + X43

# Regressions
eta5 ~ eta1 + eta2
eta6 ~ eta3 + eta4

# residual covariances
eta5 ~~ eta6
'




model2 <- SemFC$new(data=Y_2, relation_matrix = C, mode=mode, scale=F, bias=F)
model2$fit_svd()
model2$fit_ml(initialisation_svd = TRUE)


#### lavan #######
fit.sem.ml <- sem(sem.model, data=X_2)


fit <- sem(sem.model , data = X_2 , optim.gradient = "numerical")


fit <- sem(sem.model, data = X_2, verbose = TRUE)
summary (fit , standardized = T, fit.measures = T)

estimate = parameterEstimates(fit.sem.ml, standardized = TRUE)




#### csem #####
library(cSEM)
# Results for PLSc
fit_PLSc <- csem(.model = sem.model, .data = X)
summarize(fit_PLSc)
assess(fit_PLSc)
