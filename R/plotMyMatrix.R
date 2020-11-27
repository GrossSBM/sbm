#' Plot an adjacency or incidence Matrix
#'
#' @param Mat  : a matrix representing the network
#' @param dimLabels : a list of length 2 specifying the types of nodes in row and col  (functional group) (Default is \code{NULL})
#' @param clustering : a list of length 2 specifying a clustering on row and col
#' @param plotOptions : a list providing options. See details below.
#' @details The list of parameters \code{plotOptions} for the matrix plot is
#'  \itemize{
#'  \item{"legend": }{Boolean. Set TRUE if you want to see the legend. Default value is FALSE}
#'  \item{"legend.title": }{Boolean. Set TRUE if you want to print the title of the legend. Default value is FALSE}
#'  \item{"legend.position": }{Position of the legend. Possible values are 'bottom', 'top','left,'right'. Default value is 'bottom'}
#'  \item{"rowNames": }{Set true if the rownames must be plotted. Default value is FALSE}
#'  \item{"colNames": }{Set true if the colNames must be plotted. Default value is FALSE}
#'  \item{"line.color": }{Chain of character. The color of the lines to separate groups if a clustering is provided. Default value is red}
#'  \item{"line.width": }{Numeric. Width  of the lines to separate groups. Default value is NULL, automatically chosen}
#'  \item{"title": }{Chain of character. Title of the plot. Default value is NULL}
#'  }
#' @return a ggplot object corresponding to the plot
#' @export
#'
#' @examples
#' M <- matrix(sample(c(0,1),900,replace=TRUE),30,30)
#' plotMyMatrix(M, plotOptions= list(legend = FALSE))
#' M2 <- matrix( rpois(800,10),40,20)
#' plotMyMatrix(M2, dimLabels=list(row= 'reader',col = 'book'), plotOptions = list(legend = TRUE))
#'

plotMyMatrix = function(Mat, dimLabels = list(row = 'row', col = 'col'), clustering = NULL, plotOptions = NULL){

  if (is.atomic(dimLabels)){
    if (length(dimLabels) == 0){dimLabels = list(row = 'row', col = 'col')}
    if (length(dimLabels) == 1){dimLabels = list(row = dimLabels, col = dimLabels)}
    if (length(dimLabels) == 2){dimLabels = list(row = dimLabels[1], col = dimLabels[2])}
  }
  if (is.null(names(dimLabels))){names(dimLabels) = c('row','col')}

  g <- plotMatrix(Mat, dimLabels = dimLabels, clustering,plotOptions)
  g
}
