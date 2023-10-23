#' This function calculates the lambda for desired sum of absolute value
#' for the loadings
#' @param argu input vector
#' @param sumabs value of sum of absolute values for loadings desired
#' @keywords Supervised Princpal Component Analysis
#' @export
#' @examples BinarySearch(argu=u,sumabs=sqrt(10))

BinarySearch <- function(argu,sumabs){
  if(sgmeth2::l2n(argu)==0 || sum(abs(argu/sgmeth2::l2n(argu)))<=sumabs) return(0)
  lam1 <- 0
  lam2 <- max(abs(argu))-1e-5
  iter <- 1
  while(iter < 200){
    su <- sgmeth2::soft(argu,(lam1+lam2)/2)
    if(sum(abs(su/sgmeth2::l2n(su)))<sumabs){
      lam2 <- (lam1+lam2)/2
    } else {
      lam1 <- (lam1+lam2)/2
    }
    if((lam2-lam1)<1e-6) return((lam1+lam2)/2)
    iter <- iter+1
  }
  warning("Didn't quite converge")
  return((lam1+lam2)/2)
}