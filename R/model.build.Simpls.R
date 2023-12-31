#' cross validation function for sparse supervised PCA
#' @param Y vector of response 
#' @param X predictor matrix
#' @param V pca vector
#' @param kernel flag for kernel type
#' @keywords Supervised Princpal Component Analysis
#' @export
#' @examples model.build.Simpls(Y=Y,X=X,V=V,kernel="linear")

model.build.Simpls <- function(Y,X,V,kernel=c("linear","delta")){
  #Z <- Rfast::mat.mult(X,V)
  Z <- X%*%V
  df <- data.frame(y=Y,z=Z)
  if(kernel=="linear") mod <- lm(y~.,data=df)
  else{
    df$y <- factor(df$y)
    mod <- nnet::multinom(y~.,data=df,trace=F)
  }
  return(mod)
}
