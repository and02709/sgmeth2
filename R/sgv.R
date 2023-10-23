#' Calculates shrinkage for sparse group vector v
#' @param lambda group shrinkage parameter
#' @param v vector to undergo shrinkage
#' @param alpha l1 shrinkage penalty term
#' @keywords Supervised Princpal Component Analysis
#' @export
#' @examples soft(x=v,d=lambda)

sgv <- function(lambda, v, alpha){
  d <- lambda*alpha/2
  sgmeth2::l2n(sgmeth2::soft(v,d))^2 - lambda^2*(1-alpha)^2*length(v)
}