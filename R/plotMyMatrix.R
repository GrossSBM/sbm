#' Plot an adjacency or incidence Matrix
#'
#' @param Mat  : a matrix representing the network
#' @param rowLabel character : type of nodes in rows (functional group) (Default is \code{NULL})
#' @param colLabel character :  type of nodes in columns (functional group) (Default is \code{NULL})
#'
#' @return a ggplot object corresponding to the plot
#' @export
#'
#' @examples
#' M <- matrix(sample(90,c(0,1),replace=TRUE),10,9)
#' plotMyMatrix(M)
plotMyMatrix = function(Mat, rowLabel = NULL, colLabel = NULL){
  g <- plotMatrix(Mat, rowFG = rowLabel, colFG = colLabel, clustering = NULL)
  g
}
