% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotMyMatrix.R
\name{plotMyMatrix}
\alias{plotMyMatrix}
\title{Plot an adjacency or incidence Matrix}
\usage{
plotMyMatrix(
  Mat,
  dimLabels = c(row = "row", col = "col"),
  clustering = NULL,
  plotOptions = NULL
)
}
\arguments{
\item{Mat}{: a matrix representing the network}

\item{dimLabels}{: a vector of length 1 or 2 specifying the types of nodes in row and col  (functional group) (Default is \code{NULL})}

\item{clustering}{: a list of length 2 specifying a clustering on row and col}

\item{plotOptions}{: a list providing options. See details below.}
}
\value{
a ggplot object corresponding to the plot
}
\description{
Plot an adjacency or incidence Matrix
}
\details{
The list of parameters \code{plotOptions} for the matrix plot is
\itemize{
\item "legend":  Boolean. Set TRUE if you want to see the legend. Default value is FALSE
\item "legend.title":  Boolean. Set TRUE if you want to print the title of the legend. Default value is FALSE
\item "legend.position":  Position of the legend. Possible values are 'bottom', 'top','left,'right'. Default value is 'bottom'
\item "rowNames":  Set true if the rownames must be plotted. Default value is FALSE
\item "colNames":  Set true if the colNames must be plotted. Default value is FALSE
\item "line.color":  Chain of character. The color of the lines to separate groups if a clustering is provided. Default value is red
\item "line.width":  Numeric. Width  of the lines to separate groups. Default value is NULL, automatically chosen
\item "title":  Chain of character. Title of the plot. Default value is NULL
}
}
\examples{
M <- matrix(sample(c(0,1),900,replace=TRUE),30,30)
plotMyMatrix(M, dimLabels = c('individulals'), plotOptions= list(legend = FALSE))
M2 <- matrix( rpois(800,10),40,20)
plotMyMatrix(M2, dimLabels = c(row = 'reader',col = 'book'), plotOptions = list(legend = TRUE))

}
