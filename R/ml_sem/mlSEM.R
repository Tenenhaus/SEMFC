

library(Rsolnp)

source('R/ml_sem/F1.R')
source('R/ml_sem/h_constraints.R')

mlSEM <- function (init, block_sizes, mode, S, lengths_parameter,which_exo_endo){

  # number of formative blocks
  r <- sum(mode == "formative")

  if(r !=0){

  result <- solnp(pars = init,
                  fun=F1, eqfun=heq1,
                  eqB = rep(0,r), S = S, block_sizes=block_sizes, mode=mode, lengths_parameter = lengths_parameter,
                  which_exo_endo = which_exo_endo,
                  control = list(trace = 0, tol = 1e-8))
  }
  else{
    result <- solnp(pars = init,
                    fun=F1, S = S, block_sizes=block_sizes, mode=mode, lengths_parameter = lengths_parameter,
                    which_exo_endo = which_exo_endo,
                    control = list(trace = 0, tol = 1e-16))

  }

  return(result)



}