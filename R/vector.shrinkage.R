#' This function shrinks each group of the v vector by the deisred lambda penalty term
#' to achieve the desired number of nonzero groups in vector v.
#' @param Xtu.list list of estimated v vector groups with index for each loading
#' @param Xtu.vec This contains a vector of l2 norms of Xtu
#' @param Xk.col.vec This contains the number of predictors for each group
#' @param lambda group shrinkate penalty
#' @keywords Supervised Princpal Component Analysis
#' @export
#' @examples vector.shrinkage(Xtu.list = temp$Xtu.list, Xtu.vec=temp$Xtu.vec, Xk.col.vec=temp$Xk.col.vec, lambda=lambda)

vector.shrinkage <- function(Xtu.list, Xtu.vec, Xk.col.vec, lambda){
  num.group <- length(Xk.col.vec)
  new.list  <- list()
  for(k in 1:num.group){
    index <- Xtu.list[[k]][,1]
    new.vk <- sgmeth2::gsoft(pk=Xk.col.vec[k], l2=Xtu.vec[k],lambda=lambda)*Xtu.list[[k]][,2]
    new.vk.df <- data.frame(index,new.vk)
    new.list <- c(new.list,list(new.vk.df))
  }
  return(new.list)
}