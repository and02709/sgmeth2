#' cross validation function for sparse supervised PCA
#' @param X nxp predictor matrix
#' @param Y nx1 response vector
#' @param npc number of desired principal components
#' @param n.folds number of folds to perform cross validation
#' @param sparsity.type specifies mechanism for vector shrinkage
#' @param nonzero.loadings specifies the desired number of nonzero loadings
#' @param sumabsv shrinkage using sum of absolute value of loadings
#' @param kernel specification of the response kernel
#' @param parallel flag for parallel process to parApply
#' @param n.cores number of cores to be used in parallel process
#' @param niter number of iterations in the SMD algorithm
#' @param trace display algorithm progress
#' @param part.balance flag for whether folds process should balance factors
#' @param mc.method flag for whether the parallel method should use mclapply
#' @keywords Supervised Princpal Component Analysis
#' @export
#' @examples cv.SSPCA <- function(X=X, Y=Y ,npc=2, n.folds=5, sparsity.type="sumabs", nonzero.loadings=NULL, sumabsv=sqrt(10),  kernel="linear", parallel=F, n.cores=NULL, niter=200, trace=F, part.balance=F)

cv.SSPCA <- function(X, Y ,npc, n.folds, sparsity.type=c("loadings", "sumabs"), 
                     nonzero.loadings=NULL, sumabsv=NULL, 
                     kernel=c("linear", "delta"), parallel=F, n.cores=NULL, 
                     niter=50, trace=F, part.balance=T, mc.method=T){
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  
  n <- nrow(X)
  p <- ncol(X)
  
  colnames(X) <-paste0("x",c(1:p))
  colnames(Y) <- "y"
  
  if(nrow(Y)!=n) stop("number of observations in predictors and response must match")
  
  if(sparsity.type!="loadings" && sparsity.type!="sumabs") stop("Please select a valid sparsity type")
  
  if(sparsity.type=="loadings"){
    if(min(nonzero.loadings) < 1) stop("minimum number of nonzero loadings must be greater than or equal to one")
    if(max(nonzero.loadings) > p) stop("maximum number of nonzero loadings must be less than or equal to the number of predictors")
    sp.arg <- nonzero.loadings
  }
  if(sparsity.type=="sumabs"){
    if(min(sumabsv) < 1) stop("minimum number of nonzero loadings must be greater than or equal to one")
    if(max(sumabsv) > sqrt(p)) stop("maximum number of nonzero loadings must be less than or equal to the square root of the number of predictors")
    sp.arg <- sumabsv
  }
  
  if(kernel!="linear" && kernel!="delta") stop("Please select a valid kernel")
  
  df <- data.frame(Y,X)
  n.sp <- length(sp.arg)
  
  if(kernel=="delta" && part.balance){
    df.partition <- groupdata2::fold(data=df,k=n.folds,cat_col = "y")
  } else{
    df.partition <- groupdata2::fold(data=df,k=n.folds)
  }
  fold.arg <- c(1:n.folds)
  param.grid <- expand.grid(fold.arg,sp.arg)
  colnames(param.grid) <- c("fold.arg","sp.arg")
  
  if(parallel){ 
    if(mc.method){
      if(is.null(n.cores)) n.cores <- parallel::detectCores() / 2
      param.grid.l <- as.list(data.frame(t(param.grid)))
      metric.vec <- unlist(parallel::mclapply(param.grid.l,cv.partition.SSPCA,
                                              df.partition=df.partition,npc=npc,
                                              n.folds=n.folds,
                                              sparsity.type=sparsity.type,
                                              nonzero.loadings=NULL,
                                              sumabsv=NULL,kernel=kernel,
                                              niter=niter,trace=trace,
                                              mc.cores=n.cores))
    } else{
      if(is.null(n.cores)) n.cores <- parallel::detectCores() / 2
      clust <- parallel::makeCluster(n.cores)
      metric.vec <- parallel::parApply(cl=clust,X=as.matrix(param.grid),1,
                                       cv.partition.SSPCA,
                                       df.partition=df.partition,npc=npc,
                                       n.folds=n.folds,sparsity.type=sparsity.type,
                                       nonzero.loadings=NULL,
                                       sumabsv=NULL,kernel=kernel,
                                       niter=niter,trace=trace)
    }
  } else{
    metric.vec <- apply(X=as.matrix(param.grid),1,cv.partition.SSPCA,
                        df.partition=df.partition,npc=npc,
                        n.folds=n.folds,sparsity.type=sparsity.type,
                        nonzero.loadings=NULL,
                        sumabsv=NULL,kernel=kernel,niter=niter,trace=trace)
  }
  
  
  param.grid <- cbind(param.grid,metric.vec)
  metric.matrix <- sgmeth2::mat.fill(param.grid=param.grid,n.sp=n.sp,n.folds=n.folds)
  
  cv.metric <- apply(metric.matrix,1,mean)
  cv.df <- data.frame(sparsity=sp.arg,cv.metric=cv.metric)
  best.metric <- min(cv.df$cv.metric)
  best.sparse.param <- cv.df$sparsity[which(cv.df$cv.metric==best.metric)]
  
  return(list(cv.matrix=metric.matrix, cv.metrics=cv.df, best.metric=best.metric,best.sparse.param=best.sparse.param))
}