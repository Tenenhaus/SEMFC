#' BuildCall
#'
#' @param model fitted rgccac object 
#' @param display_numbers boolean if display_numbers = TRUE, 
#' display the cor_eta_y between latents and manifests nodes 
#' and the betas/gammas between the latents nodes
#' @param graph_options overlap = FALSE
#' @param node_options_LV  string. The shape and the color of the latents nodes
#' @param node_options_MV  string. The shape and the color of the manifests nodes
#' @param edge_options string. The shape and the color of the edge

#'
#' @return
#' a string which will be used as a basis for the creation of the graph
#' @export
#'
#' @examples
#' BuildCall(model = fit.rgccac, 
#' display_numbers = TRUE,
#' graph_options = list(overlap = FALSE),
#' node_options_LV = list(shape = "oval", color = "black"), 
#' node_options_MV = list(shape = "box", color = "black"),
#' edge_options = list(penwidth = 0.5, color = "black"))
BuildCall=function (model = model, 
                    display_numbers = TRUE,
                    graph_options = list(overlap = FALSE),
                    node_options_LV = list(shape = "oval", color = "black"), 
                    node_options_MV = list(shape = "box", color = "black"),
                    edge_options = list(penwidth = 0.5, color = "black"))
{
  string1 <- ""
  string1 <- paste(string1, "digraph plot {")
  string1 <- paste(string1, "\n")
  string1 <- paste(string1, "graph [", paste(paste(names(graph_options), 
                                                   graph_options, sep = " = "),
                                             collapse = ", "), 
                   "]")
  string1 <- paste(string1, "\n")
  nodes <- getNodes(model)
  string1 <- paste(string1, "node [", paste(paste(names(node_options_MV), 
                                                  node_options_MV, sep = " = "),
                                            collapse = ", "), 
                   "]")
  string1 <- paste(string1, "\n")
  string1 <- paste(string1, paste(nodes$observeds, collapse = "; "))
  string1 <- paste(string1, "\n")
  string1 <- paste(string1, "node [", paste(paste(names(node_options_LV), 
                                                  node_options_LV, sep = " = "),
                                            collapse = ", "), 
                   "]")
  string1 <- paste(string1, "\n")
  string1 <- paste(string1, paste(nodes$latents, collapse = "; "))
  string1 <- paste(string1, "\n")
  string1 <- paste(string1, "edge[", paste(paste(names(edge_options), 
                                                 edge_options, sep = " = "), 
                                           collapse = ", "), 
                   "]")
  string1 <- paste(string1, "\n")
  ifelse(display_numbers, 
         string1 <- paste(string1, model_rgcca(model, labels = TRUE)), 
         string1 <- paste(string1, model_rgcca(model, labels = FALSE))
  )
  string1 <- paste(string1, "}", sep = "\n")
  string1
}


