% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cv.gSPCA.R
\name{cv.gSPCA}
\alias{cv.gSPCA}
\title{cross validation function for groupwise supervised PCA}
\usage{
cv.gSPCA(
  X,
  Y,
  npc,
  n.folds = 5,
  groups = NULL,
  nonzero.groups = NULL,
  kernel = c("linear", "delta"),
  parallel = F,
  n.cores = NULL,
  niter = 50,
  trace = F,
  part.balance = T,
  ind.names = NULL,
  mc.method = T
)
}
\arguments{
\item{X}{nxp predictor matrix}

\item{Y}{nx1 response vector}

\item{npc}{number of desired principal components}

\item{n.folds}{number of folds to perform cross validation}

\item{groups}{vector of groups to which each predictor belongs}

\item{nonzero.groups}{number of desired nonzero groups}

\item{kernel}{specification of the response kernel}

\item{parallel}{flag for parallel process to parApply}

\item{n.cores}{number of cores to be used in parallel process}

\item{niter}{number of iterations in the SMD algorithm}

\item{trace}{display algorithm progress}

\item{part.balance}{flag for whether folds process should balance factors}

\item{ind.names}{vector of each observation label}

\item{mc.method}{flag for whether the parallel method should use mclapply}
}
\description{
cross validation function for groupwise supervised PCA
}
\examples{
cv.gSPCA(X=xtrain, Y=ytrain, npc=1, n.folds=5, groups=group.list, nonzero.groups=c(1:20), kernel="linear", parallel=F, n.cores=NULL, niter=50, trace=T, part.balance=F,ind.names=NULL)
}
\keyword{Analysis}
\keyword{Component}
\keyword{Princpal}
\keyword{Supervised}
