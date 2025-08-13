source('R/SEMFC/sem_f_c.R')


print('###### modele all reflective #############')
source('data/data_generated_reflective.R')


######  True SIGMA ########
print('DATA empirical sigma = TRUE')
#### our code #######

model <- SemFC$new(data=Y, relation_matrix = C, mode=mode, scale=F, bias=F)
model$fit_svd()
model$fit_ml(initialisation_svd = TRUE)
# model$ml_infer()



#### lavan #######

fit.sem.ml <- sem(sem.model, data=do.call(cbind, Y), estimator = "ML", likelihood="wishart")
a = fit.sem.ml@optim
f = fitMeasures(fit.sem.ml, "fmin")

estimate = parameterEstimates(fit.sem.ml, standardized = TRUE)




#### comparaison #######

lambda_th = c(l1,l2,l3,l4,l5,l6)
std_all_ml = mapply(function(x,y) sqrt(1-(x/y)), unlist(model$ml_parameters$residual_variance),diag(model$cov_S))
std_all_svd = mapply(function(x,y) sqrt(1-(x/y)), unlist(model$svd_parameters$residual_variance),diag(model$cov_S))
std_all_lavaan = estimate[1:18,11]

lambda_comparaison = cbind(estimate[1:18,1:3],lambda_th, std_all_svd, std_all_ml, std_all_lavaan)
print('lambda')
print(lambda_comparaison)


g =  unlist(model$ml_parameters$gamma)
b = unlist(model$ml_parameters$beta)
bg_ml = c(g[1,1], g[1,2], b[1,2],g[2,3],g[2,4],b[2,1])

g_svd =  unlist(model$svd_parameters$gamma)
b_svd = unlist(model$svd_parameters$beta)
bg_svd = c(g_svd[1,1], g_svd[1,2], b_svd[1,2],g_svd[2,3],g_svd[2,4],b_svd[2,1])

bg_th = c(GAMMA[1,1], GAMMA[1,2], BETA[1,2],GAMMA[2,3],GAMMA[2,4],BETA[2,1])

bg_lavaan = estimate[19:24,11]
beta_gama_comparaison = cbind(estimate[19:24,1:3], bg_th, bg_svd, bg_ml, bg_lavaan)
print('beta et gamma')
print(beta_gama_comparaison)


res_var_ml = unlist(unname(model$ml_parameters$residual_variance))
res_var_svd = unlist(unname(model$svd_parameters$residual_variance))
res_var_lavaan = estimate[26:43,11]
res_var_th = c(1-l1^2, 1-l2^2, 1-l3^2, 1-l4^2, 1-l5^2, 1-l6^2)

res_var_comparaison = cbind(estimate[25:42,1:3], res_var_th, res_var_svd, res_var_ml, res_var_lavaan)
print('residual variance')
print(res_var_comparaison)

print('F1')
f1_ml_true = model$ml_parameters$F
f1_svd_true = model$svd_parameters$F

P_exo_lavaan = lavInspect(fit.sem.ml, what = 'cor.lv')[1:4,1:4]
P_endo_lavaan = lavInspect(fit.sem.ml, what = 'cor.lv')[5:6,5:6]
g_lavaan = c(bg_lavaan[[1]], bg_lavaan[[2]], bg_lavaan[[4]], bg_lavaan[[5]])
b_lavaan = c(bg_lavaan[[3]], bg_lavaan[[6]])

param_lavaan = c(std_all_lavaan, P_exo_lavaan[upper.tri(P_exo_lavaan)], g_lavaan, b_lavaan, P_endo_lavaan[upper.tri(P_endo_lavaan)],
                 res_var_lavaan)

f1_lavaan = F1(param_lavaan, model$cov_S, model$block_sizes, model$mode, model$lengths_theta, model$which_exo_endo)

f1_true = cbind(f1_svd_true, f1_ml_true, f1_lavaan )
print(f1_true)




##### FALSE SIGMA
print('DATA empirical sigma = FALSE')

#### our code #######

model <- SemFC$new(data=Y_2, relation_matrix = C, mode=mode, scale=F, bias=F)
model$fit_svd()
model$fit_ml(initialisation_svd = TRUE)
# model$ml_infer()
# true_param_with_S
# model$summary()


#### lavan #######

fit.sem.ml <- sem(sem.model, data=do.call(cbind, Y_2), estimator = "ML", likelihood="wishart")
estimate = parameterEstimates(fit.sem.ml, standardized = TRUE)


#### comparaison #######

lambda_th = c(l1,l2,l3,l4,l5,l6)
std_all_ml = mapply(function(x,y) sqrt(1-(x/y)), unlist(model$ml_parameters$residual_variance),diag(model$cov_S))
std_all_svd = mapply(function(x,y) sqrt(1-(x/y)), unlist(model$svd_parameters$residual_variance),diag(model$cov_S))
std_all_lavaan = estimate[1:18,11]

lambda_comparaison_false = cbind(estimate[1:18,1:3],lambda_th, std_all_svd, std_all_ml, std_all_lavaan)
print('lambda')
print(lambda_comparaison_false)


g =  unlist(model$ml_parameters$gamma)
b = unlist(model$ml_parameters$beta)
bg_ml = c(g[1,1], g[1,2], b[1,2],g[2,3],g[2,4],b[2,1])

g_svd =  unlist(model$svd_parameters$gamma)
b_svd = unlist(model$svd_parameters$beta)
bg_svd = c(g_svd[1,1], g_svd[1,2], b_svd[1,2],g_svd[2,3],g_svd[2,4],b_svd[2,1])

bg_th = c(GAMMA[1,1], GAMMA[1,2], BETA[1,2],GAMMA[2,3],GAMMA[2,4],BETA[2,1])

bg_lavaan = estimate[19:24,11]
beta_gama_comparaison_false = cbind(estimate[19:24,1:3], bg_th, bg_svd, bg_ml, bg_lavaan)
print('beta et gamma')
print(beta_gama_comparaison_false)



res_var_ml = unlist(unname(model$ml_parameters$residual_variance))
res_var_svd = unlist(unname(model$svd_parameters$residual_variance))
res_var_lavaan = estimate[25:42,11]
res_var_th = c(1-l1^2, 1-l2^2, 1-l3^2, 1-l4^2, 1-l5^2, 1-l6^2)

res_var_comparaison_false = cbind(estimate[25:42,1:3], res_var_th, res_var_svd, res_var_ml, res_var_lavaan)
print('residual variance')
print(res_var_comparaison_false)


print('F1')

f1_ml_false = model$ml_parameters$F
f1_svd_false = model$svd_parameters$F

P_exo_lavaan = lavInspect(fit.sem.ml, what = 'cor.lv')[1:4,1:4]
P_endo_lavaan = lavInspect(fit.sem.ml, what = 'cor.lv')[5:6,5:6]
g_lavaan = c(bg_lavaan[[1]], bg_lavaan[[2]], bg_lavaan[[4]], bg_lavaan[[5]])
b_lavaan = c(bg_lavaan[[3]], bg_lavaan[[6]])

param_lavaan = c(std_all_lavaan, P_exo_lavaan[upper.tri(P_exo_lavaan)], g_lavaan, b_lavaan, P_endo_lavaan[upper.tri(P_endo_lavaan)],
                 res_var_lavaan)

f1_lavaan = F1(param_lavaan, model$cov_S, model$block_sizes, model$mode, model$lengths_theta, model$which_exo_endo)




f1_false = cbind(f1_svd_false,f1_ml_false, f1_lavaan)
print(f1_false)


print("######################   MIXED MODEL #####################################")


source('data/data_generated_mixed.R')

######  True SIGMA ########
print('DATA empirical sigma = TRUE')
#### our code #######

model <- SemFC$new(data=Y, relation_matrix = C, mode=mode, scale=F, bias=F)
model$fit_svd()
model$fit_ml(initialisation_svd = TRUE)


#fit.sem.ml <- sem(sem.model, data=do.call(cbind, Y), estimator = "ML", likelihood="wishart")

##### FALSE SIGMA
print('DATA empirical sigma = FALSE')

#### our code #######

model2 <- SemFC$new(data=Y_2, relation_matrix = C, mode=mode, scale=F, bias=F)
model2$fit_svd()
model2$fit_ml(initialisation_svd = TRUE)





lambda_th = c(l1,l2,l3,l4,l5,l6)
std_all_ml = unlist(model$ml_parameters$lambda)
std_all_svd =  unlist(model$svd_parameters$lambda)
std_all_ml2 = unlist(model2$ml_parameters$lambda)
std_all_svd2 = unlist(model2$svd_parameters$lambda)


lambda_comparaison3 = cbind(estimate[1:18,1:3],lambda_th, std_all_svd, std_all_ml, std_all_svd2, std_all_ml2)
print('lambda')
print(lambda_comparaison3)


g =  unlist(model$ml_parameters$gamma)
b = unlist(model$ml_parameters$beta)
bg_ml = c(g[1,1], g[1,2], b[1,2],g[2,3],g[2,4],b[2,1])

g_svd =  unlist(model$svd_parameters$gamma)
b_svd = unlist(model$svd_parameters$beta)
bg_svd = c(g_svd[1,1], g_svd[1,2], b_svd[1,2],g_svd[2,3],g_svd[2,4],b_svd[2,1])

bg_th = c(GAMMA[1,1], GAMMA[1,2], BETA[1,2],GAMMA[2,3],GAMMA[2,4],BETA[2,1])



g2 =  unlist(model2$ml_parameters$gamma)
b2 = unlist(model2$ml_parameters$beta)
bg_ml2 = c(g2[1,1], g2[1,2], b2[1,2],g2[2,3],g2[2,4],b2[2,1])

g_svd2 =  unlist(model2$svd_parameters$gamma)
b_svd2 = unlist(model2$svd_parameters$beta)
bg_svd2 = c(g_svd2[1,1], g_svd2[1,2], b_svd2[1,2],g_svd2[2,3],g_svd2[2,4],b_svd2[2,1])

bg_th = c(GAMMA[1,1], GAMMA[1,2], BETA[1,2],GAMMA[2,3],GAMMA[2,4],BETA[2,1])




beta_gama_comparaison3 = cbind(estimate[19:24,1:3], bg_th, bg_svd, bg_ml, bg_svd2, bg_ml2)
print('beta et gamma')
print(beta_gama_comparaison3)


print('F1')
f1_ml_true = model$ml_parameters$F
f1_svd_true = model$svd_parameters$F
f1_ml_false = model2$ml_parameters$F
f1_svd_false = model2$svd_parameters$F


f1_mixed = cbind(f1_svd_true, f1_ml_true, f1_svd_false,f1_ml_false)
print(f1_mixed)


print("####################### ECSI #################################")

source('data/data_ecsi.R')
#### our code #######
print('run our package')
model <- SemFC$new(data=A, relation_matrix = C_ecsi, mode=mode_ecsi, scale=F, bias=F)
model$fit_svd()
model$fit_ml(initialisation_svd = TRUE)
# model$ml_infer()


#### lavan #######
print('run lavaan')
fit.sem.ml <- sem(sem.model.ecsi, data=A2, estimator = "ML", likelihood="wishart")
estimate = parameterEstimates(fit.sem.ml, standardized = TRUE)


#### comparaison #######
print('result comparaison')
std_all_ml = mapply(function(x,y) sqrt(1-(x/y)), unlist(model$ml_parameters$residual_variance),diag(model$cov_S))
std_all_svd = mapply(function(x,y) sqrt(1-(x/y)), unlist(model$svd_parameters$residual_variance),diag(model$cov_S))
std_all_lavaan = estimate[1:18,11]

lambda_comparaison_ecsi = cbind(estimate[1:18,1:3],std_all_svd, std_all_ml, std_all_lavaan)
print('lambda')
print(lambda_comparaison_ecsi)


g =  unlist(model$ml_parameters$gamma)
b = unlist(model$ml_parameters$beta)
bg_ml = c(g[1], g[2], b[2,1],g[3],b[3,1],b[3,2],b[4,3])

g_svd =  unlist(model$svd_parameters$gamma)
b_svd = unlist(model$svd_parameters$beta)
bg_svd = c(g_svd[1], g_svd[2], b_svd[2,1],g_svd[3],b_svd[3,1],b_svd[3,2],b_svd[4,3])



bg_lavaan = estimate[19:25,11]
beta_gama_comparaison_ecsi = cbind(estimate[19:25,1:3], bg_svd, bg_ml, bg_lavaan)
print('beta et gamma')
print(beta_gama_comparaison_ecsi)


res_var_ml = unlist(unname(model$ml_parameters$residual_variance))
res_var_svd = unlist(unname(model$svd_parameters$residual_variance))
res_var_lavaan = estimate[26:43,4]


res_var_comparaison_ecsi = cbind(estimate[26:43,1:3], res_var_svd, res_var_ml, res_var_lavaan)
print('residual variance')
print(res_var_comparaison_ecsi)


print('F1')
f1_ml_ecsi = model$ml_parameters$F
f1_svd_ecsi = model$svd_parameters$F

lambda_lavaan = estimate[1:18,4]
P_exo_lavaan = lavInspect(fit.sem.ml, what = 'cor.lv')[1:1,1:1]
P_endo_lavaan = lavInspect(fit.sem.ml, what = 'cor.lv')[2:5,2:5]
g_lavaan = c(bg_lavaan[[1]], bg_lavaan[[2]], bg_lavaan[[4]])
b_lavaan = c(bg_lavaan[[3]], bg_lavaan[[5]],  bg_lavaan[[6]],  bg_lavaan[[7]])

param_lavaan = c(std_all_lavaan, g_lavaan, b_lavaan, P_endo_lavaan[upper.tri(P_endo_lavaan)],
                 res_var_lavaan)

f1_lavaan = F1(param_lavaan, model$cov_S, model$block_sizes, model$mode, model$lengths_theta, model$which_exo_endo)





f1_ecsi= cbind(f1_svd_ecsi,f1_ml_ecsi, f1_lavaan)
print(f1_ecsi)