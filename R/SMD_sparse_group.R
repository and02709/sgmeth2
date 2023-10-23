#' This function peforms penalized matrix decomposition for one pair of
#' u and v vectors
#' Uses number of nonzero groups for the shrinkage mechanism as well as l1 penalty
#' @param X matrix to undergo PMD
#' @param d supplied singular value
#' @param u supplied left vector
#' @param v supplied right vector
#' @param n number of observations
#' @param p number of predictors
#' @param groups vector of groups to which each predictor belongs
#' @param nonzero.groups number of desired nonzero groups
#' @param unique.group vector of each unique group label
#' @param num.group number of unique groups
#' @param alpha l1 shrinkage penalty
#' @param niter number of iterations for the algorithm
#' @param trace displays information about progress of algorithm
#' @param npc number of desired principal components
#' @keywords Supervised Princpal Component Analysis
#' @export
#' @examples SMD_sparse_group(X=Xuse, d=temp$d, u=temp$u, v=temp$v, n=n, p=p, groups=groups, nonzero.groups=nonzero.groups, unique.group=unique.group, alpha=alpha,num.group=num.group, niter=niter, trace=trace, npc=i)

SMD_sparse_group <- function(X, d, u, v, n, p, groups, nonzero.groups, unique.group, 
                             num.group, alpha, niter, trace, npc){
  oldv <- rnorm(p,0,1)
  oldu <- rnorm(n,0,1)
  if(trace) cat("Vector ", npc, ": ")
  for(iter in 1:niter){
    if((sum(abs(oldv-v)) < 1e-7) && (sum(abs(oldu-u)) < 1e-7)) break
    oldv <- v
    oldu <- u
    if(trace) cat(iter," ", fill=F)
    # calculate vector u
    u <- X%*%v
    # normalize vector u
    u <- matrix(matrix(u/sgmeth2::l2n(u)),ncol=1)
    # Get measurements for each group of predictors
    temp <- sgmeth2::get.Xk(X=X,u=u,groups=groups,
                   unique.group=unique.group,num.group=num.group)
    # Get lambda based on the relative strength of each group
    #   we calculate the lambda that gives us the desired number of
    #   nonzero groups
    #Xu <- t(X)%*%u
    Xu <- crossprod(X,u)
    upper.lambda <- 2*max(abs(Xu))/abs(1-alpha)+1
    upper.lambda <- min(upper.lambda, 10^5)
    
    lambda <- sgmeth2::calculate.sparse.lambda(Xtu.list=temp$Xtu.list, upper.lambda=upper.lambda,
                                      alpha=alpha, nonzero.groups=nonzero.groups)
    
    bin.list <- sgmeth2::calculate.shrinkage.binary(Xtu.list=temp$Xtu.list,lambda=lambda,alpha=alpha)
    
    # Shrink each vector by the approprite amount determined by lambda
    shrunk.list <- sgmeth2::vector.sparse.shrinkage(Xtu.list=temp$Xtu.list, shrink.list=bin.list, 
                                           lambda=lambda, alpha=alpha)
    # Replace shrunk values for each group back into the vector v
    v <- sgmeth2::embed.v(shrunk.Xtu.list=shrunk.list,vold=v)
    # normalize vector v
    v <- matrix(matrix(v/sgmeth2::l2n(v)),ncol=1)
  }
  if(trace) cat("\n")
  #d <- as.numeric(t(u)%*%(X%*%v))
  d <- as.numeric(crossprod(u,X%*%v))
  return(list(d=d, u=u, v=v))
}