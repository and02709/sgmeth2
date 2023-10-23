#' This function shrinks each group of the v vector by the deisred lambda penalty term
#' to achieve the desired number of nonzero groups in vector v while achieving l1 shrinkage.
#' @param Xtu.list list of estimated v vector groups with index for each loading
#' @param shrink.list list of groups to be zeroed out
#' @param lambda group shrinkate penalty
#' @param alpha l1 penalty term
#' @keywords Supervised Princpal Component Analysis
#' @export
#' @examples vector.sparse.shrinkage(Xtu.list=temp$Xtu.list, shrink.list=bin.list, lambda=lambda, alpha=alpha)

vector.sparse.shrinkage <- function(Xtu.list, shrink.list, lambda, alpha){
  num.group <- length(Xtu.list)
  new.list <- list()
  for(k in 1:num.group){
    index <- Xtu.list[[k]][,1]
    vk <- Xtu.list[[k]][,2]
    if(shrink.list[k]==1){
      new.vk <- rep(0, length(vk))
    } else{
      new.vk <- (1/2)*(sgmeth2::soft(vk, lambda*alpha/2) - lambda*(1-alpha)*sqrt(length(vk))*sgmeth2::soft(vk, lambda*alpha/2)/sgmeth2::l2n(sgmeth2::soft(vk, lambda*alpha/2)))
    }
    
    new.vk.df <- data.frame(index,new.vk)
    new.list <- c(new.list,list(new.vk.df))
  }
  return(new.list)
}