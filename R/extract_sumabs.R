#' This extracts and returns the left and right singular vectors
#' @param X predictor matrix
#' @param n number of observations 
#' @param p number of predictors
#' @param sumabs value of sum of absolute values for loadings desired
#' @param niter number of iterations
#' @param trace whether to output progress
#' @param npc which pc is being extracted
#' @keywords Sparse partial elast squares
#' @export
#' @examples extract_sumabs(X=X,n=200,p=1000,sumabsv=sqrt(10),niter=200,trace=F,npc=2)

extract_sumabs <- function(X, n, p, sumabsv, niter, trace, npc){
  temp <- RSpectra::svds(X, k=1)
  if(sumabsv==sqrt(p)){
    return(list(d=temp$d, u=temp$u, v=temp$v))
  } else{
    return(sgmeth2::SMD.sumabs(X=X,d=temp$d,u=temp$u,v=temp$v,n=n,p=p,
                      sumabsv=sumabsv,niter=niter,trace=trace,npc=npc))
  }
}