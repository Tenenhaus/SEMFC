% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BuildCall.R
\name{BuildCall}
\alias{BuildCall}
\title{BuildCall}
\usage{
BuildCall(
  model = model,
  display_numbers = TRUE,
  graph_options = list(overlap = FALSE),
  node_options_LV = list(shape = "oval", color = "black"),
  node_options_MV = list(shape = "box", color = "black"),
  edge_options = list(penwidth = 0.5, color = "black")
)
}
\arguments{
\item{model}{fitted rgccac object}

\item{display_numbers}{boolean if display_numbers = TRUE, 
display the cor_eta_y between latents and manifests nodes 
and the betas/gammas between the latents nodes}

\item{graph_options}{overlap = FALSE}

\item{node_options_LV}{string. The shape and the color of the latents nodes}

\item{node_options_MV}{string. The shape and the color of the manifests nodes}

\item{edge_options}{string. The shape and the color of the edge}

\item{C}{matrix of 1 and 0 with the same number of row and column without 
any 1 on the diag}
}
\value{
a string which will be used as a basis for the creation of the graph
}
\description{
BuildCall
}
\examples{
BuildCall(model = fit.rgccac, 
display_numbers = TRUE,
graph_options = list(overlap = FALSE),
node_options_LV = list(shape = "oval", color = "black"), 
node_options_MV = list(shape = "box", color = "black"),
edge_options = list(penwidth = 0.5, color = "black"))
}
