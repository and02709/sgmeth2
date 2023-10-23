#' This function calculates the lambda necessary to achieve the desired number
#' of nonzero groups.
#' @param Xtu.vec This contains a vector of l2 norms of Xtu
#' @param Xk.col.vec This contains the number of predictors for each group
#' @param nonzero.groups number of desired nonzero groups
#' @keywords Supervised Princpal Component Analysis
#' @export
#' @examples calculate.lambda(Xtu.vec=temp$Xtu.vec,Xk.col.vec=temp$Xk.col.vec,nonzero.groups=nonzero.groups)

calculate.lambda <- function(Xtu.vec, Xk.col.vec, nonzero.groups){
  num.group <- length(Xk.col.vec)
  unorder.vec <- rep(0, num.group)
  for(k in 1:num.group){
    unorder.vec[k] <- 2*Xtu.vec[k]/sqrt(Xk.col.vec[k])
  }
  index.vec <- order(unorder.vec, decreasing = T)
  order.vec <- unorder.vec[index.vec]
  return(sum(order.vec[nonzero.groups]+order.vec[nonzero.groups+1])/2)
}