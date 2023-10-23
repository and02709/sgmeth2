#' This function peforms penalized matrix decomposition for one pair of
#' u and v vectors
#' Uses sum of absolute value of loadings as the shrinkage mechanism
#' @param X matrix to undergo PMD
#' @param d supplied singular value
#' @param u supplied left vector
#' @param v supplied right vector
#' @param n number of observations
#' @param p number of predictors
#' @param sumasv sum of absolute value of loadings of the v vector
#' @param niter number of iterations for the algorithm
#' @param trace displays information about progress of algorithm
#' @param npc number of desired principal components
#' @keywords Supervised Princpal Component Analysis
#' @export
#' @examples SMD.sumabs(X=X,d=D,u=U,v=V,n=200,p=1000,sumabsv=sqrt(10),niter=200,trace=F,npc=1)

SMD.sumabs <- function(X, d, u, v, n, p, sumabsv, niter, trace, npc){
  oldv <- rnorm(p,0,1)
  oldu <- rnorm(n,0,1)
  if(trace) cat("Vector ", npc, ": ")
  for(iter in 1:niter){
    if((sum(abs(oldv-v)) < 1e-7) && (sum(abs(oldu-u)) < 1e-7)) break
    oldv <- v
    oldu <- u
    if(trace) cat(iter," ", fill=F)
    # update u
    #argu <- Rfast::mat.mult(X,v)
    argu <- X%*%v
    u <- matrix(argu/sgmeth2::l2n(argu),ncol=1)
    # update v
    #argv <- Rfast::Crossprod(u,X)
    argv <- crossprod(u,X)
    #v <- matrix(argv/l2n(argv),ncol=1)
    # Find appropriate shrinkage 
    lamv <- sgmeth2::BinarySearch(argu=argv,sumabs=sumabsv)
    # soft threshold v
    sv <- sgmeth2::soft(argv,lamv)
    v <- matrix(sv/sgmeth2::l2n(sv),ncol=1)
  }
  if(trace) cat("\n")
  #d <- as.numeric(t(u)%*%(X%*%v))
  d <- as.numeric(crossprod(u,X%*%v))
  return(list(d=d, u=u, v=v))
}
