#' cross validation function for groupwise SIMPLS
#' @param X nxp predictor matrix
#' @param Y nx1 response vector
#' @param npc number of desired principal components
#' @param n.folds number of folds to perform cross validation
#' @param kernel specification of the response kernel
#' @param groups vector of groups to which each predictor belongs
#' @param nonzero.groups number of desired nonzero groups
#' @param trace display algorithm progress
#' @param niter number of iterations in the SMD algorithm
#' @param parallel flag for parallel process to parApply
#' @param n.cores number of cores to be used in parallel process
#' @param ind.names vector of each observation label
#' @param part.balance flag for whether folds process should balance factors
#' @param mc.method flag for whether the parallel method should use mclapply
#' @keywords Supervised Princpal Component Analysis
#' @export
#' @examples cv.gSIMPLS(X, Y ,npc, n.folds, groups=NULL, nonzero.groups=NULL, kernel=c("linear", "delta"), parallel=F, n.cores=NULL, niter=50, trace=F, part.balance=T,ind.names=NULL)

cv.gSIMPLS <- function(X, Y ,npc, n.folds, groups=NULL, nonzero.groups=NULL, 
                       kernel=c("linear", "delta"), parallel=F, n.cores=NULL, 
                       niter=50, trace=F, part.balance=T,mc.method=T,ind.names=NULL){
  
  # Convert X to a matrix  
  X <- as.matrix(X)
  # Convert Y to a matrix
  Y <- as.matrix(Y)
  # Find number of observations
  n <- dim(X)[1]
  # Find number of predictors
  p <- dim(X)[2]
  # Stop if number of observations doesn't match between X and Y
  if(nrow(Y)!=n) stop("number of observations in predictors and response must match")
  
  # Check names of predictors
  if(is.null(colnames(X))) colnames(X) <- paste0("X",seq(1:p))
  colnames(Y) <- "y"
  
  # label rows of X and Y
  if(is.null(ind.names)){
    rownames(X) <- c(1:n)
    rownames(Y) <- c(1:n)
    ind.names <- rownames(X)
  }
  if(length(ind.names)!=dim(X)[[1]] || length(ind.names)!= dim(Y)[[1]]){
    stop("ind.names length does not match number of observations in data")
  }
  
  if(kernel!="linear" && kernel!="delta") stop("Please select a valid kernel")
  
  if(is.null(groups)) stop("Please supply group information for each predictor in X")
  if(length(groups) != p) stop("Group labels do not match number of predictors")
  unique.groups <- unique(groups)
  num.group <- length(unique.groups)
  
  if(sum((nonzero.groups-floor(nonzero.groups))==0) != length(nonzero.groups)) stop("Must specify integer groups")
  
  if(min(nonzero.groups) < 1) stop("Must specify minimum number of nonzero groups as at least 1")
  if(max(nonzero.groups) > num.group) stop("Cannot have more nonzero groups than total number of groups")
  
  df <- data.frame(Y,X)
  num.nz.gr <- length(nonzero.groups)
  
  if(kernel=="delta" && part.balance){
    df.partition <- groupdata2::fold(data=df,k=n.folds,cat_col = "y")
  } else{
    df.partition <- groupdata2::fold(data=df,k=n.folds)
  }
  fold.arg <- c(1:n.folds)
  param.grid <- expand.grid(fold.arg,nonzero.groups)
  colnames(param.grid) <- c("fold.arg","gr.arg")
  
  if(parallel){ 
    if(mc.method){
      if(is.null(n.cores)) n.cores <- parallel::detectCores() 
      param.grid.l <- as.list(data.frame(t(param.grid)))
      metric.vec <- unlist(parallel::mclapply(param.grid.l,cv.partition.group.simpls,
                                              df.partition=df.partition,npc=npc, 
                                              n.folds=n.folds, groups=groups, 
                                              kernel=kernel, 
                                              ind.names=ind.names, 
                                              niter=niter,trace=trace,
                                              mc.cores=n.cores))
    } else{
      if(is.null(n.cores)) n.cores <- parallel::detectCores() 
      clust <- parallel::makeCluster(n.cores)
      metric.vec <- parallel::parApply(cl=clust,X=as.matrix(param.grid),1,
                                       cv.partition.group.simpls,df.partition=df.partition,
                                       npc=npc, n.folds=n.folds, groups=groups, 
                                       kernel=kernel, ind.names=ind.names, 
                                       niter=niter,trace=trace) 
    }
  } else{
    metric.vec <- apply(X=as.matrix(param.grid),1,cv.partition.group.simpls,
                        df.partition=df.partition,npc=npc, n.folds=n.folds, 
                        groups=groups, kernel=kernel, ind.names=ind.names, 
                        niter=niter,trace=trace)
  }
  
  
  param.grid <- cbind(param.grid,metric.vec)
  metric.matrix <- sgmeth2::mat.fill(param.grid=param.grid,n.sp=num.nz.gr,n.folds=n.folds)
  
  cv.metric <- apply(metric.matrix,1,mean)
  cv.df <- data.frame(nonzero.groups=nonzero.groups,cv.metric=cv.metric)
  best.metric <- min(cv.df$cv.metric)
  best.group.param <- cv.df$nonzero.groups[which(cv.df$cv.metric==best.metric)]
  
  return(list(cv.matrix=metric.matrix, cv.metrics=cv.df, best.metric=best.metric,best.group.param=best.group.param))
}
