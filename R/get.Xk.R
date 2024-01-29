#' This function calculates the v vector loadings, l2 norms for each group, and number of preditors in each group
#' @param X predictor matrix
#' @param u iterated vector corresponding to observations
#' @param groups vector of groups to which each predictor belongs
#' @param unique.group list of unique group labesl
#' @param num.group number of unique groups
#' @keywords Supervised Princpal Component Analysis
#' @export
#' @examples get.Xk(X=X,u=u,groups=groups,unique.group=unique.group,num.group=num.group)

get.Xk <- function(X,u,groups,unique.group,num.group){
  Xtu.vec <- rep(0,num.group)
  Xk.col.vec <- rep(0,num.group)
  Xtu.list <- list()
  for(k in 1:num.group){
    index <- which(groups==unique.group[k])
    X.k <- X[,index]
    Xtu <- t(X.k)%*%u
    Xtu.df <- data.frame(index,Xtu)
    
    Xtu.list <- c(Xtu.list,list(Xtu.df))
    if(is.null(ncol(X.k))){
      Xk.col <- 1
      } else{
      Xk.col <- ncol(X.k)
      }
    Xtu.vec[k] <- sgmeth2::l2n(Xtu)
    Xk.col.vec[k] <- Xk.col
  }
  return(list(Xtu.list=Xtu.list,Xtu.vec=Xtu.vec,Xk.col.vec=Xk.col.vec))
}
