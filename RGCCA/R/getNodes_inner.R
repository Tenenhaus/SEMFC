#'GetNodes
#' Give the names of the latents variables
#'
#' @param fit a fitted rgccac object
#'
#' @return
#' list of names of the latents variables
#' @export
#'
#' @examples
getNodes_inner = function(fit){
  
  latent_nodes = names(fit$d)
  list( latents = latent_nodes)
}