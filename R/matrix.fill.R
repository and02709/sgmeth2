#' cross validation function for sparse supervised PCA
#' @param param.grid
#' @param n.args number of grid search parameter options
#' @param n.folds number of folds
#' @keywords Supervised Princpal Component Analysis
#' @export
#' @examples matrix.fill(param.grid=param.grid,n.args=n.args,n.folds=n.folds)


matrix.fill <- function(param.grid,n.args,n.folds){
  #metric.matrix <- matrix(NA,nrow=n.lams,ncol = n.folds)
  back.set <- n.folds-1
  metric.matrix <- sapply(1:n.args,function(x) param.grid[(x*n.folds-back.set):(x*n.folds),4])
  metric.matrix <- t(metric.matrix)
  return(metric.matrix)
}