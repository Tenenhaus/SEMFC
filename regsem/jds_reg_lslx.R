library(lslx)
library(regsem)
library(lavaan)


source('R/SEMFC/sem_f_c.R')
source('regsem/utils_function_sparse.R')
source('regsem/regsem_model.R')

Ns <- c(300, 500, 1000)





table_nul_FP <- matrix(NA, nrow = length(Ns), ncol = 2)


table_fort_FN <- matrix(NA, nrow = length(Ns), ncol = 2)

table_faible_FN <- matrix(NA, nrow = length(Ns), ncol = 2)




colnames(table_nul_FP) <- colnames(table_fort_FN) <- colnames(table_faible_FN) <-
      c('regsem', 'lslx')

rownames(table_nul_FP) <- rownames(table_fort_FN) <- rownames(table_faible_FN)  <- paste0('N=', Ns)



table_time <- matrix(NA, nrow = length(Ns), ncol = 2)
colnames(table_time) <- c('regsem', 'lslx')
rownames(table_time) <- paste0('N=', Ns)


# range_low <- 131:160
# range_strong <- 1:30
# range_null <- 31:130

range_low <- 41:60
range_strong <- 1:20
range_null <- 21:40



len_block <- 60
models <- regssem_model(len_block)
sem.model.regsem <- models$sem.model
sem.model.lslx <- models$sem.model.lslx




for (N in Ns){
  print(N)
  source('data/data_generated_reflective.R')
  X <- X_2
  Y <- Y_2


  lslx_fa <- lslx$new(model = sem.model.lslx, data = X)

  start <- proc.time()

  lslx_fa$fit(
    penalty_method = "mcp",
    lambda_grid = exp(seq(log(0.001), log(1), length.out = 10)),
    delta_grid = c(1, 2, 3, 5, 10, Inf)
  )
  end <- proc.time()
  time_lslx <- end - start


  start <- proc.time()
  lav <- sem(sem.model.regsem, X)
  regsem.out <- cv_regsem(lav, type="lasso", pars_pen = 1:(len_block-1),n.lambda=23,jump=.05)
  end <- proc.time()
  time_regsem <- end - start



  raw_lambda_lslx = rowSums(lslx_fa$extract_coefficient_matrix(selector = "bic", block = 'y<-f')$g)
  raw_lambda_lslx = raw_lambda_lslx[colnames(X)]

  raw_lambda_regsem = regsem.out$final_pars[grep("-> X", names(regsem.out$final_pars))]
  names(raw_lambda_regsem) <- sub(".*->\\s*", "", names(raw_lambda_regsem))
  raw_lambda_regsem <- c(raw_lambda_regsem, X11 = 1, X21 = 1, X31=1, X41=1, X51=1, X61=1)[colnames(X)]






  fnr_faible_regsem = classification(raw_lambda_regsem[range_low], l1[range_low])$FNR
  fnr_fort_regsem = classification(raw_lambda_regsem[range_strong], l1[range_strong])$FNR
  fpr_nul_regsem = classification(raw_lambda_regsem[range_null], l1[range_null])$FPR

  fnr_faible_lslx = classification(raw_lambda_lslx[range_low], l1[range_low])$FNR
  fnr_fort_lslx = classification(raw_lambda_lslx[range_strong], l1[range_strong])$FNR
  fpr_nul_lslx = classification(raw_lambda_lslx[range_null], l1[range_null])$FPR




  table_fort_FN[which(Ns == N), ] <- c(fnr_fort_regsem, fnr_fort_lslx)

  table_faible_FN[which(Ns == N), ] <- c(fnr_faible_regsem, fnr_faible_lslx)

  table_nul_FP[which(Ns == N), ] <- c(fpr_nul_regsem, fpr_nul_lslx)



  table_time[which(Ns == N), ] <- c(time_regsem[3], time_lslx[3])



  print(table_nul_FP)
  print(table_fort_FN)
  print(table_faible_FN)
  print(table_time)

}
