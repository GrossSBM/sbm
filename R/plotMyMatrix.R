#' Plot an adjacency or incidence Matrix
#'
#' @param Mat  : a matrix representing the network
#' @param dimLabels : a list of length 2 specifying the types of nodes in row and col  (functional group) (Default is \code{NULL})
#' @param clustering : a list of length 2 specifying a clustering on row and col
#' @param plotOptions: a list. plotOptions$legend= FALSE if legend must not be given
#' @return a ggplot object corresponding to the plot
#' @export
#'
#' @examples
#' M <- matrix(sample(c(0,1),900,replace=TRUE),30,30)
#' plotMyMatrix(M)
#' M2 <- matrix(rpois(800,10),40,20)
#' plotMyMatrix(M2,dimLabels=list(row= 'reader',col = 'book'))
#'

plotMyMatrix = function(Mat, dimLabels = list(row = NULL, col = NULL),clustering = NULL,plotOptions=NULL){

  if (is.atomic(dimLabels)){
    if (length(dimLabels) == 0){dimLabels = list(row = NULL, col = NULL)}
    if (length(dimLabels) == 1){dimLabels = list(row = dimLabels, col =dimLabels )}
    if (length(dimLabels) == 2){dimLabels = list(row = dimLabels[1], col = dimLabels[2])}
  }
  if (is.null(names(dimLabels))){names(dimLabels) = as.list(dimLabels)}

  g <- plotMatrix(Mat, dimLabels = dimLabels, clustering,plotOptions)
  g
}
