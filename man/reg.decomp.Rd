% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reg.decomp.R
\name{reg.decomp}
\alias{reg.decomp}
\title{This function decomposes the learned matrix}
\usage{
reg.decomp(X, L, npc, n, p)
}
\arguments{
\item{X}{nxp matrix of predictors}

\item{L}{nxn response kernel}

\item{npc}{number of desired principal components}

\item{n}{an integer storing the number response observations}

\item{p}{an integer storing the number of predictors}
}
\description{
This function decomposes the learned matrix
}
\examples{
reg.decomp(X=X,L=L,npc=2,n=200,p=1000)
}
\keyword{Analysis}
\keyword{Component}
\keyword{Princpal}
\keyword{Supervised}
