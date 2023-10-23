#' This extracts and returns the left and right singular vectors
#' @param X predictor matrix
#' @param n number of observations 
#' @param p number of predictors
#' @param load value for number of desired nonzero loadings
#' @param niter number of iterations
#' @param trace whether to output progress
#' @param npc which pc is being extracted
#' @keywords Sparse partial elast squares
#' @export
#' @examples extract_load(X=X,n=200,p=1000,nonzero.loadings=10,niter=200,trace=F,npc=2)

extract_load <- function(X, n, p, nonzero.loadings, niter, trace, npc){
  temp <- RSpectra::svds(X, k=1)
  if(nonzero.loadings==p){
    return(list(d=temp$d, u=temp$u, v=temp$v))
  } else{
    return(sgmeth2::SMD.load(X=X,d=temp$d,u=temp$u,v=temp$v,n=n,p=p,
                    nonzero.loadings=nonzero.loadings,
                    niter=niter,trace=trace,npc=npc))
  }
}
