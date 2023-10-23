#' cross validation function for groupwise SIMPLS
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
#' @examples cv.partition.group.simpls(arg.group, df.partition, npc, n.folds, groups, kernel=c("linear", "delta"), ind.names, niter,trace)

cv.partition.group.simpls <- function(arg.group, df.partition, npc, n.folds, groups, kernel=c("linear", "delta"), ind.names, niter,trace){
  test.index <- arg.group[[1]]
  nonzero.groups <- arg.group[[2]]
  
  xtrain <- df.partition |> dplyr::filter(.folds!=test.index) |> dplyr::ungroup() |> dplyr::select(-c(y,.folds))
  xtest <- df.partition |> dplyr::filter(.folds==test.index) |> dplyr::ungroup() |> dplyr::select(-c(y,.folds))
  ytrain <- df.partition |> dplyr::filter(.folds!=test.index) |> dplyr::ungroup() |> dplyr::select(y)
  ytest <- df.partition |> dplyr::filter(.folds==test.index) |> dplyr::ungroup() |> dplyr::select(y)
  
  xmeans <- colMeans(xtrain)
  xtrain <- t(apply(xtrain, 1, function(x) x-xmeans))
  xtest <- t(apply(xtest, 1, function(x) x-xmeans))
  
  gr.simpls.out <- sgmeth2::gSIMPLS(X=xtrain,Y=ytrain,npc=npc, kernel=kernel, groups=groups, nonzero.groups=nonzero.groups, ind.names=rownames(xtrain), niter=niter, trace=trace)
  
  V <- gr.simpls.out
  
  mod <- sgmeth2::model.build.Simpls(Y=as.matrix(ytrain),X=as.matrix(xtrain),V=V,kernel=kernel)
  metric <- sgmeth2::model.metric.Simpls(Y=as.matrix(ytest),X=as.matrix(xtest),V=V,mod=mod,kernel=kernel)
  return(metric)
}