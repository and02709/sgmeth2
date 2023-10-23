#' This function calculates the lambda necessary to achieve the desired number
#' of nonzero groups while incorporating an l1 penalty on the vector
#' @param Xtu.list list containing each group for vector v and index for each loading
#' @param upper.lambda maximum lambda allowed for uniroot function
#' @param alpha l1 sparseness penalty
#' @param nonzero.groups desired number of nonzero groups remaining in vector v
#' @keywords Supervised Princpal Component Analysis
#' @export
#' @examples calculate.sparse.lambda(Xtu.list=Xtu.list, upper.lambda=upper.lambda, alpha=alpha, nonzero.groups=nonzero.groups)

calculate.sparse.lambda <- function(Xtu.list, upper.lambda, alpha, nonzero.groups){
  num.group <- length(Xtu.list)
  unorder.vec <- rep(0, num.group)
  for(k in 1:num.group){
    vk <- Xtu.list[[k]][,2]
    min.lambda.k <- uniroot(sgmeth2::sgv, lower=0, upper=upper.lambda, v=vk, alpha=alpha)
    unorder.vec[k] <- min.lambda.k$root
  }
  order.vec <- unorder.vec[order(unorder.vec, decreasing = T)]
  return(sum(order.vec[nonzero.groups]+order.vec[nonzero.groups+1])/2)
}