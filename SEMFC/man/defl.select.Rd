% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defl.select.R
\name{defl.select}
\alias{defl.select}
\title{deflation function}
\usage{
defl.select(yy, rr, nncomp, nn, nbloc)
}
\arguments{
\item{yy}{A matrix that contains the SGCCA block components of each block: \eqn{y_{1h}, \ldots, y_{Jh}}}

\item{rr}{A list that contains the residual matrices \eqn{X_{1h}, \ldots, X_{Jh}}}

\item{nncomp}{A \eqn{1 \times J} vector that contains the number of components to compute for each block.}

\item{nn}{A \eqn{1 \times J} vector that contains the numbers of already computed components for each block}

\item{nbloc}{Number of blocks.}
}
\value{
\item{resdefl}{A list of \eqn{J} elements that contains \eqn{X_{1,h+1}, \ldots, X_{J,h+1}}.}

\item{pdefl}{A list of \eqn{J} elements that contains \eqn{p_{1h}, \ldots, p_{Jh}}.}
}
\description{
The function defl.select() computes residual matrices \eqn{X_{1,h+1}, \ldots, X_{J,h+1}}. These 
residual matrices are determined according to the following formula: \eqn{X_{j,h+1} = X_{jh} - y_{jh} p_{jh}^t}.
}
