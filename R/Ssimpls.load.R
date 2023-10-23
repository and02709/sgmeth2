#' This function performs sparse partial least squares using number of nonzero loadings.
#' @param X nxp predictor matrix
#' @param Y nx1 response vector
#' @param npc number of desired principal components
#' @param nonzero.loadings desired number of nonzero loadings
#' @param kernel specification of the response kernel
#' @param niter number of iterations in the SMD algorithm
#' @param trace display algorithm progress
#' @keywords Sparse partial least squares
#' @export
#' @examples Ssimpls.load(X=X,Y=Y,npc=2,nonzero.loadings=10,kernel="linear",niter=200,trace=F)

Ssimpls.load <- function(X, Y, npc, 
                         nonzero.loadings=NULL, 
                         kernel=c("linear", "delta"), niter=50, trace=F){
  
  # Find number of observations
  n <- dim(X)[1]
  # Find number of predictors
  p <- dim(X)[2]
  
  # Check nonzero.loadings
  if(is.null(nonzero.loadings)) nonzero.loadings <- p
  if(nonzero.loadings < 1) stop("Please provide a number of nonzero loadings between 1 and p")
  if(nonzero.loadings > p) stop("Please provide a number of nonzero loadings between 1 and p")
  
  # Check names of predictors
  if(is.null(colnames(X))) colnames(X) <- paste0("X",seq(1:p))
  
  # Generate response kernel
  L <- sgmeth2::respkernel(Y=Y, n=n, kernel=kernel)
  # Initialize vector containing vectors
  v <- matrix(0, nrow=p, ncol=npc)
  Xi <- X
  
  for(i in 1:npc){
    Psi <- sgmeth2::get.psi(X=Xi, L=L, n=n, p=p)
    udv <- sgmeth2::extract_load(X=Psi, n=n, p=p, nonzero.loadings=nonzero.loadings, 
                            niter=niter, trace=trace, npc=i)
    v[,i] <- udv$v
    if(i==npc) break
    Xi <- sgmeth2::deflate(X=X,V=v[,1:i],n=n,p=p)
  }
  
  return(v)
}