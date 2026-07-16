devtools::load_all()





source('inst/model/model_reflective.R')
Y <- Y_2
X <- X_2

set.seed(27)
modelsvd <- SemFC$new(data=Y, relation_matrix = C, mode=mode, estimator = "svd")

modelsvd$fit(infer = F)





init <- modelsvd$get_estimate('theta')
model <- modelsvd$get_model()
S <- cov(X)
r = 4
tol = 1e-08

f <- function (x) { return(F1(x, S, model)) }
grad <- function(x){ return(grad_F1(x, S, model)) }
h_eq <- function(x){ return(heq(x, S, model)) }
grad_h_eq <- function(x){ return(grad_heq(x, S, model)) }

system.time({



result <- solnp(pars = init,
                    fun=F1, S = S, model = model,
                    control = list(trace = 0, tol = tol))
})






system.time({

res<- nloptr(
      x0=init,
      eval_f=f,
      eval_grad_f=grad,
      opts = list("algorithm"="NLOPT_LD_LBFGS",
                  'xtol_rel' = tol)
    )
})




system.time({
sol <- csolnp(pars = init, fn =f , gr = grad, lower = -10*abs(init), upper = 10*abs(init),
              control = list(trace = 0, tol = tol), use_r_version = T)
})





round(result$pars - res$solution, 3)
round(result$pars - sol$pars, 3)



