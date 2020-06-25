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
#' M <- matrix(sample(c(0,1),900,replace=TRUE),30,30)
#' plotMyMatrix(M)
#' M2 <- matrix(rpois(800,10),40,20)
#' plotMyMatrix(M2,rowLabel='ind',colLabel = 'book')
#'

plotMyMatrix = function(Mat, rowLabel = NULL, colLabel = NULL){
  g <- plotMatrix(Mat, rowFG = rowLabel, colFG = colLabel, clustering = NULL)
  g
}
