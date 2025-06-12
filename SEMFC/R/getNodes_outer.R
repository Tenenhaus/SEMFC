#' Outer GetNodes
#'
#' Give the names of the latents and manifests variables
#'
#' @param fit  a fitted rgccac object
#' @param num string or integer. Display the outer graph of the num block or
#' all the blocks if num  = "all"
#'
#' @return
#' list of names of the latents and manifests variables
#' @export
#'
#' @examples
getNodes_outer = function(fit, num){
  if(num == "all"){
    observed_nodes = as.vector(unlist(sapply(fit$lambda, rownames)))
    latent_nodes = names(fit$d)
    list(observeds = observed_nodes, latents = latent_nodes)
  }else{
    observed_nodes = as.vector(unlist(sapply(fit$lambda[num], rownames)))
    latent_nodes = names(fit$d[num])
    list(observeds = observed_nodes, latents = latent_nodes)
  }
}