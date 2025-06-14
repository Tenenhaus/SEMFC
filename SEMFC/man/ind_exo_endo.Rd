% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ind_exo_endo.R
\name{ind_exo_endo}
\alias{ind_exo_endo}
\title{ind_exo_endo}
\usage{
ind_exo_endo(C)
}
\arguments{
\item{C}{matrix filled with 1 and 0 with the same 
number of columns and rows and without any 1 on the diagonal}
}
\value{
list of H, J, index of the endogenous variables and 
index of the exogenous variables of the C matrix
}
\description{
ind_exo_endo
}
\examples{
C = matrix(c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,1,
0,0,1,1,1,0),6,6)
ind_exo_endo(C)
}
