% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Ssimpls.load.R
\name{Ssimpls.load}
\alias{Ssimpls.load}
\title{This function performs sparse partial least squares using number of nonzero loadings.}
\usage{
Ssimpls.load(
  X,
  Y,
  npc,
  nonzero.loadings = NULL,
  kernel = c("linear", "delta"),
  niter = 50,
  trace = F
)
}
\arguments{
\item{X}{nxp predictor matrix}

\item{Y}{nx1 response vector}

\item{npc}{number of desired principal components}

\item{nonzero.loadings}{desired number of nonzero loadings}

\item{kernel}{specification of the response kernel}

\item{niter}{number of iterations in the SMD algorithm}

\item{trace}{display algorithm progress}
}
\description{
This function performs sparse partial least squares using number of nonzero loadings.
}
\examples{
Ssimpls.load(X=X,Y=Y,npc=2,nonzero.loadings=10,kernel="linear",niter=200,trace=F)
}
\keyword{Sparse}
\keyword{least}
\keyword{partial}
\keyword{squares}
