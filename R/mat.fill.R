#' cross validation function for sparse supervised PCA
#' @param param.grid data frame containing parameters
#' @param n.sp number of sparse parameter options
#' @param n.folds number of folds
#' @keywords Supervised Princpal Component Analysis
#' @export
#' @examples mat.fill(param.grid=param.grid,n.sp=100,n.folds=5)


mat.fill <- function(param.grid,n.sp,n.folds){
  metric.matrix <- matrix(NA, nrow=n.sp,ncol=n.folds)
  for(i in 1:n.sp){
    beg.index <- i*n.folds - (n.folds-1)
    last.index <- i*n.folds
    metric.matrix[i,] <- param.grid[beg.index:last.index,3]
  }
  return(metric.matrix)
}
