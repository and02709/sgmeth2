#' This function replaces the oldings of the old v vector with the new values
#' @param shrunk.Xtu.list list of shrunk loadings for each group with corresponding index
#' @param vold old v vector to be replaced with new loadings
#' @keywords Supervised Princpal Component Analysis
#' @export
#' @examples embed.v(shrunk.Xtu.list=shrunk.list,vold=v)

embed.v <- function(shrunk.Xtu.list, vold){
  num.group <- length(shrunk.Xtu.list)
  for(k in 1:num.group){
    temp.df <- shrunk.Xtu.list[[k]]
    vold <- replace(vold,temp.df$index,temp.df$new.vk)
  }
  return(vold)
}