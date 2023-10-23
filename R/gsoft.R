#' This function performs group soft thresholding on a vector.
#' @param pk kth group of vector
#' @param l2 l2 norm
#' @param lambda shrinkage penalty term
#' @keywords Supervised Princpal Component Analysis
#' @export
#' @examples soft(x=v,d=lambda)

gsoft <- function(pk,l2,lambda){
  return(max(0, (1-(lambda/2)*sqrt(pk)/l2)))
}