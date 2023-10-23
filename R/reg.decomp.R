#' This function decomposes the learned matrix
#' 
#' @param X nxp matrix of predictors
#' @param L nxn response kernel
#' @param npc number of desired principal components
#' @param n an integer storing the number response observations
#' @param p an integer storing the number of predictors 
#' @keywords Supervised Princpal Component Analysis
#' @export
#' @examples reg.decomp(X=X,L=L,npc=2,n=200,p=1000)

reg.decomp <- function(X, L, npc, n, p){
  H <- diag(x=1, nrow=n) - 1/n*rep(1,n)%*%t(rep(1,n))
  #M <- Rfast::Crossprod(X, Rfast::mat.mult(H, Rfast::mat.mult(L, Rfast::mat.mult(H, X))))
  M <- crossprod(X, H%*%L%*%H%*%X)
  Md <- RSpectra::eigs_sym(M, k=npc, which = "LM")
  return(list(vectors=Md$vectors, values=Md$values))
}
