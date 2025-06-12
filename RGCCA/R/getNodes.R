


#' GetNodes
#' Give the names of the latents and manifests variables
#'
#' @param fit a fitted rgccac object
#'
#' @return list of names of the latents and manifests variables
#' @export
#'
#' @examples
#' getNodes(fit.rgccac)

getNodes = function(fit){
  observed_nodes = as.vector(unlist(sapply(fit$lambda, rownames)))
  latent_nodes = names(fit$d)
  list(observeds = observed_nodes, latents = latent_nodes)
}




