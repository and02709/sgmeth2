#' This function performs sparse partial least squares using sum of absolute values of loadings.
#' @param X nxp predictor matrix
#' @param Y nx1 response vector
#' @param npc number of desired principal components
#' @param sumabsv shrinkage using sum of absolute value of loadings
#' @param kernel specification of the response kernel
#' @param niter number of iterations in the SMD algorithm
#' @param trace display algorithm progress
#' @keywords Supervised Princpal Component Analysis
#' @export
#' @examples Ssimpls.sumabs(X=X,Y=Y,npc=2,sumabsv=sqrt(10),kernel="linear",niter=200,trace=F)


Ssimpls.sumabs <- function(X, Y, npc, 
                           sumabsv=NULL, 
                           kernel=c("linear", "delta"), niter=50, trace=F){
  
  # Find number of observations
  n <- dim(X)[1]
  # Find number of predictors
  p <- dim(X)[2]
  
  # Check nonzero.loadings
  if(is.null(sumabsv)) sumabsv <- sqrt(p)
  if(sumabsv < 1) stop("Please provide a number of nonzero loadings between 1 and sqrt(p)")
  if(sumabsv > sqrt(p)) stop("Please provide a number of nonzero loadings between 1 and sqrt(p)")
  
  # Check names of predictors
  if(is.null(colnames(X))) colnames(X) <- paste0("X",seq(1:p))
  
  # Generate response kernel
  L <- sgmeth2::respkernel(Y=Y, n=n, kernel=kernel)
  # Initialize vector containing vectors
  v <- matrix(0, nrow=p, ncol=npc)
  Xi <- X
  
  for(i in 1:npc){
    Psi <- sgmeth2::get.psi(X=Xi, L=L, n=n, p=p)
    udv <- sgmeth2::extract_sumabs(X=Psi, n=n, p=p, sumabsv=sumabsv, niter=niter, 
                              trace=trace, npc=i)
    v[,i] <- udv$v
    if(i==npc) break
    Xi <- sgmeth2::deflate(X=X,V=v[,1:i],n=n,p=p)
  }
  
  return(v)
}