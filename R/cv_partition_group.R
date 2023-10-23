#' cross validation function for groupwise supervised PCA
#' @param arg.group vector for folds and sparse parameter
#' @param df.partition list of folds
#' @param npc number of desired principal components
#' @param n.folds number of folds to perform cross validation
#' @param groups vector of groups to which each predictor belongs
#' @param kernel specification of the response kernel
#' @param ind.names vector of each observation label
#' @param niter number of iterations in the SMD algorithm
#' @param trace display algorithm progress
#' @keywords Supervised Princpal Component Analysis
#' @export
#' @examples cv_partition_group(arg.group=c(test.index,sparse.lambda), df.partition=groupdata2.partition,npc=2,n.folds=5,kernel="linear",ind.names=NULL,niter=200,trace=F)

cv_partition_group <- function(arg.group, df.partition, npc, n.folds, groups, kernel=c("linear", "delta"), ind.names, niter,trace){
  test.index <- arg.group[[1]]
  nonzero.groups <- arg.group[[2]]
  
  xtrain <- df.partition |> dplyr::filter(.folds!=test.index) |> dplyr::ungroup() |> dplyr::select(-c(y,.folds))
  xtest <- df.partition |> dplyr::filter(.folds==test.index) |> dplyr::ungroup() |> dplyr::select(-c(y,.folds))
  ytrain <- df.partition |> dplyr::filter(.folds!=test.index) |> dplyr::ungroup() |> dplyr::select(y)
  ytest <- df.partition |> dplyr::filter(.folds==test.index) |> dplyr::ungroup() |> dplyr::select(y)
  
  xmeans <- colMeans(xtrain)
  xtrain <- t(apply(xtrain, 1, function(x) x-xmeans))
  xtest <- t(apply(xtest, 1, function(x) x-xmeans))
  
  gr.spca.out <- gSPCA(X=xtrain,Y=ytrain,npc=npc, kernel=kernel, groups=groups, nonzero.groups=nonzero.groups, ind.names=rownames(xtrain), niter=niter, trace=trace)
  
  V <- gr.spca.out$v
  
  mod <- model.build.SSPCA(Y=as.matrix(ytrain),X=as.matrix(xtrain),V=V,kernel=kernel)
  metric <- model.metric.SSPCA(Y=as.matrix(ytest),X=as.matrix(xtest),V=V,mod=mod,kernel=kernel)
  return(metric)
}