
library(lavaan)
library(devtools)
load_all()

print("####################### ECSI #################################")

source('inst/model/model_ecsi.R')
#### our code #######
print('run our package')
model <- SemFC$new(data=A, relation_matrix = C_ecsi, mode=mode_ecsi, scale=F, bias=F, estimator = "svd")
modelml <- SemFC$new(data=A, relation_matrix = C_ecsi, mode=mode_ecsi, scale=F, bias=F)


model$fit(infer = F)
modelml$fit(infer = F)

#### lavan #######
print('run lavaan')
fit.sem.ml <- sem(sem.model.ecsi, data=A2, estimator = "ML", likelihood="wishart" )
estimate = parameterEstimates(fit.sem.ml, standardized = TRUE)


#### comparaison #######
print('result comparaison')


std_all_ml = unlist(modelml$estimate$std_lambda)
std_all_svd = unlist(model$estimate$std_lambda)


std_all_lavaan = estimate[1:18,11]

lambda_comparaison_ecsi = cbind(estimate[1:18,1:3],std_all_svd, std_all_ml, std_all_lavaan)
print('lambda')
print(lambda_comparaison_ecsi)


g =  unlist(modelml$estimate$gamma)
b = unlist(modelml$estimate$beta)
bg_ml = c(g[1], g[2], b[2,1],g[3],b[3,1],b[3,2],b[4,3])

g_svd =  unlist(model$estimate$gamma)
b_svd = unlist(model$estimate$beta)
bg_svd = c(g_svd[1], g_svd[2], b_svd[2,1],g_svd[3],b_svd[3,1],b_svd[3,2],b_svd[4,3])



bg_lavaan = estimate[19:25,11]
beta_gama_comparaison_ecsi = cbind(estimate[19:25,1:3], bg_svd, bg_ml, bg_lavaan)
print('beta et gamma')
print(beta_gama_comparaison_ecsi)


res_var_ml = unlist(unname(modelml$estimate$residual_variance))
res_var_svd = unlist(unname(model$estimate$residual_variance))
res_var_lavaan = estimate[26:43,4]


res_var_comparaison_ecsi = cbind(estimate[26:43,1:3], res_var_svd, res_var_ml, res_var_lavaan)
print('residual variance')
print(res_var_comparaison_ecsi)


print('F1')
f1_ml_ecsi = modelml$gof$F
f1_svd_ecsi = model$gof$F

lambda_lavaan = estimate[1:18,4]
P_exo_lavaan = lavInspect(fit.sem.ml, what = 'cor.lv')[1:1,1:1]
P_endo_lavaan = lavInspect(fit.sem.ml, what = 'cor.lv')[2:5,2:5]
g_lavaan = c(bg_lavaan[[1]], bg_lavaan[[2]], bg_lavaan[[4]])
b_lavaan = c(bg_lavaan[[3]], bg_lavaan[[5]],  bg_lavaan[[6]],  bg_lavaan[[7]])

param_lavaan = c(std_all_lavaan, g_lavaan, b_lavaan,
                 res_var_lavaan)

f1_lavaan = F1(param_lavaan, model$model$cov_S, model$model)





f1_ecsi= cbind(f1_svd_ecsi,f1_ml_ecsi, f1_lavaan)
print(f1_ecsi)