library(R6)


# load_all("SEMFC")
source("SEMFC/R/svdSEM.R")
source("SEMFC/R/scale2.R")
source("SEMFC/R/correction.R")
source("SEMFC/R/cov2.R")
source("SEMFC/R/lvm.R")
source("SEMFC/R/ind_exo_endo.R")
source("SEMFC/R/d_LS.R")
source("SEMFC/R/improper.R")
source("SEMFC/R/svdSEM_infer.R")
source("SEMFC/R/scaleDataSet.R")
source("SEMFC/R/svdSEM_gof.R")
source("SEMFC/R/model_sem.R")
source("SEMFC/R/mlSEM_infer.R")

#functions
source("functions/data_simulation.R")
source("functions/F1.R")
source("functions/h_theta.R")
source("functions/s_implied.R")
source("functions/h_constraints.R")
source("functions/parameters_svd.R")
source("functions/mlSEM.R")




library(Matrix)
library(knitr)
library(pheatmap)

library(MASS)
library(lavaan)
library(Rsolnp)




SEM_F_C <- R6Class(
  "SEM_F_C",
  public = list(
    # Attributs
    data = NULL,
    relation_matrix = NULL,
    which_exo_endo = NULL,
    scale = FALSE,
    mode = NULL,
    cov_S = NULL,
    bias = FALSE,
    svd_result = NULL,
    n_blocks = NULL,
    n_row = NULL,
    block_sizes = NULL,
    lengths_theta = NULL,
    S_composites = NULL,
    boot_svd = NULL,
    ml_infer_estimate = NULL,
    SD = NULL,
    VCOV = NULL,
    svd_parameters = list(),
    ml_parameters = list(),







    # Méthode d'initialisation
    initialize = function(data, relation_matrix, mode,scale, bias) {
      self$data <- data
      self$relation_matrix <- relation_matrix
      self$scale <- ifelse(is.null(scale), FALSE, scale)
      self$bias <- ifelse(is.null(bias), FALSE, bias)
      self$mode <- mode

      parameter_model <- get_parameter_model_sem(data, mode)
      self$n_blocks <- parameter_model$n_blocks
      self$n_row <- parameter_model$n_row
      self$block_sizes <- parameter_model$block_sizes
      self$cov_S <- parameter_model$S
      self$S_composites <- parameter_model$S_diag_composites

      which_exo_endo <- ind_exo_endo(relation_matrix)
      self$which_exo_endo <-which_exo_endo

      self$lengths_theta <- get_lengths_theta(relation_matrix, self$block_sizes, mode)




    },



    # Méthode fit utilisant la technique SVD
    fit_svd = function() {
      svd_result <- svdSEM(self$data,
                           self$relation_matrix,
                           self$scale,
                           self$mode,
                           self$bias)

      self$svd_parameters <- svd_result
      theta_svd <- parameters_svd(lambda = svd_result$lambda,
                                  P_EXO = svd_result$P_EXO,
                                  G = svd_result$gamma,
                                  B = svd_result$beta,
                                  P_ENDO = svd_result$P_ENDO,
                                  residual_variance = svd_result$residual_variance,
                                  S_composites = self$S_composites,
                                  mode = self$mode)

      self$svd_parameters$theta <- theta_svd
      self$svd_parameters$F <- F1_bis(theta_svd, self$cov_S, self$block_sizes, self$mode, self$lengths_theta, self$which_exo_endo)


    },
    svd_infer = function(B = 1000, verbose = TRUE){
      if (is.null(self$svd_parameters)) {
        self$fit_svd()
      }
      boot_out <- svdSEM_infer(self$svd_parameters, B = 1000, verbose = TRUE)
      self$boot_svd <- boot_out
    },


    fit_ml = function(initialisation_svd = TRUE) {

      block_sizes <- self$block_sizes
      C <- self$relation_matrix
      mode <- self$mode

      # Initialisation par SVD si demandé
      if (initialisation_svd) {
        if (is.null(self$svd_parameters$theta)) {
          self$fit_svd() }

        initial_params <- self$svd_parameters$theta
      } else {
        initial_params <- NULL
      }

      ml_sol <- mlSEM(initial_params, block_sizes, mode, self$cov_S, self$lengths_theta, self$which_exo_endo)
      theta_ml <- ml_sol$pars
      self$ml_parameters <- s_implied_bis( x = theta_ml, block_sizes = block_sizes, mode =mode,
                                           lengths_parameter = self$lengths_theta, which_exo_endo = self$which_exo_endo, jac = F)
      self$ml_parameters$theta <- theta_ml
      self$ml_parameters$F <- F1_bis(theta_ml, self$cov_S, self$block_sizes, self$mode, self$lengths_theta, self$which_exo_endo)


    },



    ml_infer = function(){
      theta_ml <- self$ml_parameters$theta
      block_sizes <- self$block_sizes
      mode <- self$mode
      S <- self$cov_S
      N <- self$n_row



      ml_infer_estimate <- mlSEM_infer(theta_ml, S, block_sizes, mode, self$lengths_theta, N, self$ml_parameters,self$which_exo_endo)

      self$ml_infer_estimate <- ml_infer_estimate
    },



    summary = function(){
      cat("\nCoefficients:\n")
      cat("lambda:\n")
      printCoefmat(model$ml_infer_estimate$estimate$lambda, P.values = TRUE, has.Pvalue = TRUE)
      cat("beta:\n")
      printCoefmat(model$ml_infer_estimate$estimate$beta, P.values = TRUE, has.Pvalue = TRUE)
      cat("gamma:\n")
      printCoefmat(model$ml_infer_estimate$estimate$gamma, P.values = TRUE, has.Pvalue = TRUE)


    }



  )
)


library(readxl)
ECSI = as.data.frame(read_excel("mobil.xls"))/10
A = list(CUSTOMER_E = ECSI[, c("CUEX1", "CUEX2", "CUEX3")],
         PERC_QUAL  = ECSI[, c("PERQ1", "PERQ2", "PERQ3", "PERQ4",
                               "PERQ5", "PERQ6", "PERQ7")],
         PERC_VALUE = ECSI[, c("PERV1", "PERV2")],
         CUSTOMER_S = ECSI[, c("CUSA1", "CUSA2", "CUSA3")],
         CUSTOMER_L = ECSI[, c("CUSL1", "CUSL2", "CUSL3")])

C <- matrix(c(0, 1, 1, 1, 0,
              0, 0, 1, 1, 0,
              0, 0, 0, 1, 0,
              0, 0, 0, 0, 1,
              0, 0, 0, 0, 0), 5, 5, byrow = TRUE)


colnames(C) <- rownames(C) <- names(A)

mode <- c(rep("reflective", 5))


model <- SEM_F_C$new(data=A, relation_matrix = C, mode=mode, scale=F, bias=F)
model$fit_svd()
x = model$svd_parameters$theta
res =s_implied_bis(x, model$block_sizes, model$mode, model$lengths_theta, model$which_exo_endo, jac=F)


# model$svd_infer(B=1000, verbose=FALSE)
model$fit_ml(initialisation_svd = TRUE)
# model$ml_infer()





A2 = data.frame(Reduce("cbind", A))

sem.model <-  '
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

fit.sem.ml <- sem(sem.model, data=A2, estimator = "ML")


estimate = parameterEstimates(fit.sem.ml, standardized = TRUE)




print('end')






set.seed(20091979)
N <- 300
library(mvtnorm)
X <- mvrnorm(N, rep(0, 18), SIGMA, empirical = FALSE)
colnames(X) <- paste("X", rep(1:6, each = 3), rep(1:3, 6), sep ="")


Y <- list(X1 = X[, 1:3], X2 = X[, 4:6], X3 = X[, 7:9],
         X4 = X[, 10:12], X5 = X[, 13:15], X6 = X[, 16:18])



C <- matrix(c(0, 0, 0, 0, 1, 0,
             0, 0, 0, 0, 1, 0,
             0, 0, 0, 0, 0, 1,
             0, 0, 0, 0, 0, 1,
             0, 0, 0, 0, 0, 1,
             0, 0, 0, 0, 1, 0), 6, 6, byrow = TRUE)

colnames(C) <- rownames(C) <- names(Y)

mode <- c(rep("formative", 4), rep("reflective", 2))

model <- SEM_F_C$new(data=Y, relation_matrix = C, mode=mode, scale=F, bias=F)
model$fit_svd()
# model$svd_infer(B=1000, verbose=FALSE)
model$fit_ml(initialisation_svd = TRUE)
model$ml_infer()
print(model$SD)



print("end")

