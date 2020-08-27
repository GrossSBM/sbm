#' Plot an adjacency or incidence Matrix
#'
#' @param Mat  : a matrix representing the network
#' @param dimLabels : a list of length 2 specifying the types of nodes in row and col  (functional group) (Default is \code{NULL})
#' @return a ggplot object corresponding to the plot
#' @export
#'
#' @examples
#' M <- matrix(sample(c(0,1),900,replace=TRUE),30,30)
#' plotMyMatrix(M)
#' M2 <- matrix(rpois(800,10),40,20)
#' plotMyMatrix(M2,dimLabels=list(row= 'reader',col = 'book'))
#'

plotMyMatrix = function(Mat, dimLabels = list(row = NULL, col = NULL)){
  g <- plotMatrix(Mat, dimLabels = dimLabels, clustering = NULL)
  g
}
