library(lavaan)
library(lslx)
library(regsem)
library(lessSEM)
library(PMA)




source('data/data_generated_reflective.R')
source('R/SEMFC/sem_f_c.R')


pen_vars <- paste0("X1", 2:20, collapse = " + ")
pen_line <- paste0("pen() * eta1 =~ ", pen_vars)



sem.model <-paste0('
# latent variable definitions
eta1 =~ X11+',pen_vars,'
eta2 =~ X21+X22+X23
eta3 =~ X31+X32+X33
eta4 =~ X41+X42+X43
eta5 =~ X51+X52+X53
eta6 =~ X61+X62+X63

# Regressions
eta5 ~ eta1 + eta2 + eta6
eta6 ~ eta3 + eta4 + eta5

# residual covariances
eta5 ~~ eta6')


sem.model.lslx <-  paste0('
# latent variable definitions
eta1 =~ X11
eta2 =~ X21+X22+X23
eta3 =~ X31+X32+X33
eta4 =~ X41+X42+X43
eta5 =~ X51+X52+X53
eta6 =~ X61+ X62+X63

# penalisation
',pen_line,'

# Regressions
eta5 ~ eta1 + eta2 + eta6
eta6 ~ eta3 + eta4 + eta5

#variance

eta1 ~~ 1*eta1
eta2 ~~ 1*eta2
eta3 ~~ 1*eta3
eta4 ~~ 1*eta4
eta5 ~~ 1*eta5
eta6 ~~ 1*eta6




# residual covariances
eta5 ~~ eta6


')


sem.model.lsem <-paste0('
# latent variable definitions
eta1 =~ X11+',pen_vars,'
eta2 =~ X21+X22+X23
eta3 =~ X31+X32+X33
eta4 =~ X41+X42+X43
eta5 =~ X51+X52+X53
eta6 =~ X61+X62+X63

# Regressions
eta5 ~ eta1 + eta2 + eta6
eta6 ~ eta3 + eta4 + eta5

# residual covariances
eta5 ~~ eta6

')

X = X_2
Y = Y_2

model <- SemFC$new(data=Y, relation_matrix = C, mode=mode, scale=F, bias=F, pen=F,
                   len_seq = NULL,
                   nfold=NULL,
                   niter=NULL,
                   sparse_val=NULL)
model$fit_svd()
model$fit_ml(initialisation_svd = TRUE)

lav <- sem(sem.model, X)
estimate = parameterEstimates(lav, standardized = TRUE)


lslx_fa <- lslx$new(model = sem.model.lslx, data = X)

start <- proc.time()

lslx_fa$fit(
  penalty_method = "mcp",
  lambda_grid = exp(seq(log(0.001), log(1), length.out = 10)),
  delta_grid = c(1, 1.5, 2, 3, 5, 10, Inf)
)
end <- proc.time()
elapsed <- end - start
print('time lslx:')
print(elapsed)


start <- proc.time()
regsem.out <- cv_regsem(lav, type="lasso", pars_pen = 1:19,n.lambda=23,jump=.05)
end <- proc.time()
elapsed <- end - start
print('time regsem:')
print(elapsed)

regularised <- paste0("eta1=~X1", 2:20)
lavaanModel <- sem(sem.model.lsem, data = X)

start <- proc.time()
lsem <- lasso(
  lavaanModel = lavaanModel,
  regularized =regularised,
  nLambdas = 50)

end <- proc.time()
elapsed <- end - start
print('time lsem:')
print(elapsed)



#formate lavaan stand lambda
lambda_lavaan  = estimate[estimate$op == "=~",c(1:3,11) ]

#formate lslx stand lambda
raw_lambda_lslx = rowSums(lslx_fa$extract_coefficient_matrix(selector = "bic", block = 'y<-f')$g)

# common_names <- intersect(names(raw_lambda_lslx), names(var_reg_lslx))

lambda_lslx = mapply(function(x, y) x/sqrt(y), raw_lambda_lslx[colnames(X)], diag(cov(X)))



#formate regsem stand lambda
raw_lambda_regsem = regsem.out$final_pars[grep("-> X", names(regsem.out$final_pars))]
names(raw_lambda_regsem) <- sub(".*->\\s*", "", names(raw_lambda_regsem))
raw_lambda_regsem <- c(raw_lambda_regsem, X11 = 1, X21 = 1, X31=1, X41=1, X51=1, X61=1)[colnames(X)]

raw_variances_regsem <- regsem.out$final_pars[grep("~~ e", names(regsem.out$final_pars))]
is_var <- sapply(strsplit(names(raw_variances_regsem), " ~~ "), function(x) x[1] == x[2])
raw_variances_regsem <- raw_variances_regsem[is_var]
raw_variances_regsem <- unlist(mapply(rep, raw_variances_regsem, lengths(lambda), SIMPLIFY = FALSE))

lambda_regsem <- mapply(function(x, y, z) x*sqrt(y)/sqrt(z), raw_lambda_regsem, raw_variances_regsem, diag(cov(X)))


#formate lsem stand lambda

estimate_lsem = estimates(lsem, criterion = "BIC")
raw_lambda_lsem = estimate_lsem[,grep("=~X", colnames(estimate_lsem))]
names(raw_lambda_lsem) <- sub(".*=~\\s*", "", names(raw_lambda_lsem))
raw_lambda_lsem <- c(raw_lambda_lsem, X11 = 1, X21 = 1, X31=1, X41=1, X51=1, X61=1)[colnames(X)]

raw_variances_lsem <- estimate_lsem[, grep("~~e", colnames(estimate_lsem))]
is_var_lsem <- sapply(strsplit(names(raw_variances_lsem), "~~"), function(x) x[1] == x[2])
raw_variances_lsem <- raw_variances_lsem[is_var_lsem]
raw_variances_lsem <- unlist(mapply(rep, raw_variances_lsem, lengths(lambda), SIMPLIFY = FALSE))


lambda_lsem <- mapply(function(x, y, z) x*sqrt(y)/sqrt(z), raw_lambda_lsem, raw_variances_lsem, diag(cov(X)))



lambda_th = c(l1,l2,l3,l4,l5,l6)
std_all_ml = mapply(function(x,y) sqrt(1-(x/y)), unlist(model$ml_parameters$residual_variance),diag(model$cov_S))
std_all_svd = mapply(function(x,y) sqrt(1-(x/y)), unlist(model$svd_parameters$residual_variance),diag(model$cov_S))


modelpen <- SemFC$new(data=Y, relation_matrix = C, mode=mode, scale=F, bias=F, pen='pmd',
                      len_seq = NULL,
                      nfold=NULL,
                      niter=NULL,
                      sparse_val=3.1)
modelpen$fit_svd()

pensvd = mapply(function(x,y) sqrt(1-(x/y)), unlist(modelpen$svd_parameters$residual_variance),diag(modelpen$cov_S))

table_reg_empirical = cbind(lambda_th, lambda_lavaan, std_all_ml, std_all_svd, pensvd, lambda_lslx, lambda_regsem, lambda_lsem)
# table_lambda_empirical = cbind(lambda_th, lambda_lavaan, std_all_ml, std_all_svd, lambda_lslx, lambda_regsem, lambda_lsem)


