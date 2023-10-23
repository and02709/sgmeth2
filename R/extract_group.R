#' This extracts and returns the left and right singular vectors using shrinkage on the predictor vector by groups
#' @param X predictor matrix
#' @param n number of observations 
#' @param p number of predictors
#' @param groups vector containing groups for which each predictor belongs
#' @param nonzero.groups number of desired nonzero groups
#' @param unique.group vector containing unique group listings
#' @param num.group number of unique groups
#' @param niter number of iterations
#' @param trace whether to output progress
#' @param npc which pc is being extracted
#' @keywords Sparse partial elast squares
#' @export
#' @examples extract_group(X, n, p, groups, nonzero.groups, unique.group, num.group, niter, trace, npc)

extract_group <- function(X, n, p, groups, nonzero.groups, unique.group, num.group, niter, trace, npc){
  temp <- RSpectra::svds(X, k=1)
  return(sgmeth2::SMD_group(X=X,d=temp$d,u=temp$u,v=temp$v,n=n,p=p,
                    groups=groups, nonzero.groups=nonzero.groups, unique.group=unique.group,
                    num.group=num.group, niter=niter,trace=trace,npc=npc))
}
