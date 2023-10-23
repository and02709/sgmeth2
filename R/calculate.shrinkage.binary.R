#' This function calculates which groups should be completely zeroed
#' @param Xtu.list list containing each group for vector v and index for each loading
#' @param lambda minimum lambda that achieves desired number of nonzero groups
#' @param alpha l1 sparseness penalty
#' @keywords Supervised Princpal Component Analysis
#' @export
#' @examples calculate.shrinkage.binary(Xtu.list=Xtu.list,lambda=lambda,alpha=alpha)

calculate.shrinkage.binary <- function(Xtu.list,lambda,alpha){
  num.group <- length(Xtu.list)
  d1 <- (1/2)*lambda*alpha
  d2 <- lambda*(1-alpha)
  bin.shrink <- rep(0, num.group)
  for(k in 1:num.group){
    vk <- Xtu.list[[k]][,2]
    bin.shrink[k] <- (sgmeth2::l2n(sgmeth2::soft(vk,d1)) <= d2*sqrt(length(vk)))
  }
  return(bin.shrink)
}