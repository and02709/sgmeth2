#' This function peforms penalized matrix decomposition on a supplied matrix
#' and uses nonzero groups as the shrinkage mechanism
#' @param X matrix to undergo PMD
#' @param npc number of desired principal components
#' @param n number of observations
#' @param p number of predictors
#' @param groups vector of groups to which each predictor belongs
#' @param nonzero.groups number of desired nonzero groups
#' @param unique.group vector of each unique group label
#' @param num.group number of unique groups
#' @param niter number of iterations for the algorithm
#' @param trace displays information about progress of algorithm
#' @keywords Supervised Princpal Component Analysis
#' @export
#' @examples PMD_group(X=Psi,npc=npc,n=n,p=p,groups=groups,nonzero.groups=nonzero.groups,unique.group=unique.group,num.group=num.group,niter=niter,trace=trace)

PMD_group <- function(X,npc,n,p,groups,nonzero.groups,unique.group,
                      num.group,niter,trace){
  # Initialize the list of singular values and orthogonal matrices
  d <- rep(0, npc)
  u <- matrix(0, nrow=n, ncol=npc)
  v <- matrix(0, nrow=p, ncol=npc)
  
  # Initialize X matrix to be used
  Xuse <- X
  
  for(i in 1:npc){
    temp <- RSpectra::svds(Xuse, k=1)
    if(i==1){
      temp <- sgmeth2::SMD_group(X=Xuse, d=temp$d, u=temp$u, v=temp$v, n=n, p=p, 
                        groups=groups, nonzero.groups=nonzero.groups, 
                        unique.group=unique.group, 
                        num.group=num.group,niter=niter, trace=trace, npc=i)
    }
    else{
      temp <- sgmeth2::SMDorth_group(X=Xuse, us=u[,1:(i-1)], d=temp$d, u=temp$u, 
                            v=temp$v, n=n, p=p, 
                            groups=groups, nonzero.groups=nonzero.groups, 
                            unique.group=unique.group, 
                            num.group=num.group,niter=niter, trace=trace, npc=i)
    }
    d[i] <- temp$d
    u[,i] <- temp$u
    v[,i] <- temp$v
    
    if(i==npc) break
    Xuse <- Xuse - d[1]*u[,1]%*%t(v[,1])
    
  }
  
  return(list(d=d,u=u,v=v))
  
}