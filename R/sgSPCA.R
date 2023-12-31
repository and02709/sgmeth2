#' This function is the parent function for performing sparse groupwise supervised PCA
#' @param X nxp matrix of predictors
#' @param Y nx1 response vector
#' @param npc number of desired principal components
#' @param kernel designation for either linear or delta kernel
#' @param groups vector of groups to which each predictor belongs
#' @param nonzero.groups number of desired nonzero groups
#' @param ind.names vector of each observation label
#' @param alpha l1 penalty term
#' @param niter number of iterations for the algorithm
#' @param trace displays information about progress of algorithm
#' @keywords Supervised Princpal Component Analysis
#' @export
#' @examples sgSPCA(X=X,Y=Y,npc=2, kernel="linear", groups=group.list, nonzero.groups=5, ind.names=NULL, alpha=0.5, niter=200, trace=F)

sgSPCA <- function(X,Y,npc, kernel=c("linear","delta"), groups=NULL, 
                   nonzero.groups=NULL, ind.names=NULL, alpha=0,
                   niter=50, trace=F){
  X <- as.matrix(X)
  # Convert Y to a matrix
  Y <- as.matrix(Y)
  # Stop if number of observations doesn't match between X and Y
  if(dim(X)[1] != dim(Y)[1]) stop("Number of observations do not match")
  # Find number of observations
  n <- dim(X)[1]
  # Find number of predictors
  p <- dim(X)[2]
  
  # Check names of predictors
  if(is.null(colnames(X))) colnames(X) <- paste0("X",seq(1:p))
  
  # label rows of X and Y
  if(is.null(ind.names)){
    rownames(X) <- c(1:n)
    rownames(Y) <- c(1:n)
    ind.names <- rownames(Y)
  }
  if(length(ind.names)!=dim(X)[[1]] || length(ind.names)!= dim(Y)[[1]]){
    stop("ind.names length does not match number of observations in data")
  }
  
  # Generate response kernel
  L <- sgmeth2::respkernel(Y=Y, n=n, kernel=kernel)
  
  # Perform SPCA if no groups are specified and nonzero.groups is not specified
  if(is.null(groups) && is.null(nonzero.groups)){
    return(sgmeth2::SSPCA(X=X,Y=Y,npc=npc,nonzero.loadings=p, sparsity.type="loadings",
                 kernel=kernel, niter=niter, trace=trace))
  }
  
  # Perform SSPCA if no groups are specified but values of nonzero.groups is specified
  if(is.null(groups) && !is.null(nonzero.groups)){
    return(sgmeth2::SSPCA(X=X,Y=Y,npc=npc,nonzero.loadings=nonzero.groups, sparsity.type="loadings",
                 kernel=kernel, niter=niter, trace=trace))
  }
  
   if(alpha==0 && !is.null(groups) && !is.null(nonzero.groups)){
      return(sgmeth2::gSPCA(X=X,Y=Y,npc=npc, kernel=kernel, groups=groups, 
                 nonzero.groups=nonzero.groups, ind.names=ind.names, 
                 niter=niter, trace=trace) )
  }
  
  # Determine number of unique groups
  unique.group <- unique(groups)
  num.group <- length(unique(groups))
  if(!(nonzero.groups%%1==0)) stop("Must specify integer number of nonzero groups")
  if(nonzero.groups < 1) stop("Must specify at least 1 nonzero group")
  if(nonzero.groups > num.group) stop("Cannot have more nonzero groups than total number of groups")
  if(nonzero.groups == num.group){
    test.penalty <- sqrt(p)*(1-alpha)
    sparse.penalty <- test.penalty
    if(test.penalty < 1) sparse.penalty <- 1
    return(sgmeth2::SSPCA(X=X,Y=Y,npc=npc,sumabsv = sparse.penalty, sparsity.type="sumabs", kernel=kernel, niter=niter, trace=trace))
  }
  
  # if(nonzero.groups <=0 || nonzero.groups >= num.group){
  #   browser()
  #   stop("Incorrect number of nonzero groups")
  # } 
  
  # Construct the Psi matrix to perform group decomposition
  Psi <- sgmeth2::get.psi(X=X, L=L, n=n, p=p)
  
  # Group decomposition
  out <- sgmeth2::sparse.group.PMD(X=Psi,npc=npc,n=n,p=p,groups=groups,
                          nonzero.groups=nonzero.groups,alpha=alpha,
                          unique.group=unique.group,num.group=num.group,
                          niter=niter,trace=trace)
  
  # return relevant vectors
  return(out)
}
