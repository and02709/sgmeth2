#' cross validation function for sparse supervised PCA
#' @param X nxp predictor matrix
#' @param Y nx1 response vector
#' @param npc number of desired principal components
#' @param n.folds number of folds to perform cross validation
#' @param kernel specification of the response kernel
#' @param groups vector of groups to which each predictor belongs
#' @param nonzero.groups number of desired nonzero groups
#' @param alpha l1 penalty term
#' @param trace display algorithm progress
#' @param niter number of iterations in the SMD algorithm
#' @param parallel flag for parallel process to parApply
#' @param n.cores number of cores to be used in parallel process
#' @param ind.names vector of each observation label
#' @param part.balance flag for whether folds process should balance factors
#' @param mc.method flag for whether the parallel method should use mclapply
#' @keywords Supervised Princpal Component Analysis
#' @export
#' @examples cv.gSPCA(X=xtrain, Y=ytrain, npc=1, n.folds=5, groups=group.list, nonzero.groups=c(1:20), kernel="linear", parallel=F, n.cores=NULL, niter=50, trace=T, part.balance=F,ind.names=NULL)

cv.sgSPCA <- function(X,Y,npc, n.folds=5, kernel=c("linear","delta"), groups=NULL, 
                      nonzero.groups=NULL, alpha=NULL, trace=F,niter=200,
                      parallel=F,n.cores=NULL, ind.names=NULL, part.balance=F, mc.method=T){
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
  
  if(is.null(alpha)){
    alpha <- seq(from=0,to=1,by=0.05)
    num.alpha.step <- length(alpha)
    alpha[num.alpha.step] <- 0.99
  } 
  
  if(min(alpha < 0)) stop("Must have non-negative alpha values")
  if(max(alpha >= 1)) stop("Must have alpha values less than 1")
  
  df <- data.frame(Y,X)
  
  if(kernel=="delta" && part.balance){
    df.partition <- groupdata2::fold(data=df,k=n.folds,cat_col = "y")
  } else{
    df.partition <- groupdata2::fold(data=df,k=n.folds)
  }
  
  n.args <- length(nonzero.groups)*length(alpha)
  fold.arg <- c(1:n.folds)
  param.grid <- expand.grid(fold.arg,nonzero.groups,alpha)
  colnames(param.grid) <- c("fold.arg","nonzero.groups","alpha")
  
  if(parallel){ 
    if(mc.method){
      if(is.null(n.cores)) n.cores <- parallel::detectCores() 
      param.grid.l <- as.list(data.frame(t(param.grid)))
      metric.vec <- unlist(parallel::mclapply(param.grid.l,cv.partition.sparse_group,
                                              df.partition=df.partition,npc=npc,
                                              n.folds=n.folds, groups=groups, 
                                              kernel=kernel, ind.names=ind.names, 
                                              niter=niter,trace=trace,
                                              mc.cores = n.cores))
    } else{
      if(is.null(n.cores)) n.cores <- parallel::detectCores() 
      clust <- parallel::makeCluster(n.cores)
      metric.vec <- parallel::parApply(cl=clust,X=as.matrix(param.grid),1,
                                       cv.partition.sparse_group,
                                       df.partition=df.partition,npc=npc,
                                       n.folds=n.folds, groups=groups, 
                                       kernel=kernel, ind.names=ind.names, 
                                       niter=niter,trace=trace)
    }
  } else{
    metric.vec <- apply(X=as.matrix(param.grid),1,cv.partition.sparse_group,
                        df.partition=df.partition,npc=npc,
                        n.folds=n.folds, groups=groups, kernel=kernel, 
                        ind.names=ind.names, niter=niter,trace=trace)
  }
  
  param.grid <- cbind(param.grid,metric.vec)
  metric.matrix <- sgmeth2::matrix.fill(param.grid=param.grid,n.args=n.args,n.folds=n.folds)
  cv.metric <- rowMeans(metric.matrix)
  arg.list.vec <- expand.grid(nonzero.groups,alpha)
  cv.df <- data.frame(arg.list.vec,cv.metric)
  colnames(cv.df) <- c("nonzero.groups","alpha","cv.metric")
  best.metric <- min(cv.df$cv.metric)
  best.nonzero.groups <- cv.df$nonzero.groups[which(cv.df$cv.metric==best.metric)]
  best.alpha <- cv.df$alpha[which(cv.df$cv.metric==best.metric)]
  
  return(list(cv.matrix=metric.matrix, cv.metric=cv.df, nonzero.groups=best.nonzero.groups, alpha=best.alpha))
}
