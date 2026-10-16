#' Plot an alluvial plot between  clusterings
#'
#' @param listMemberships  : a list vectors containing the memberships
#' @param plotOptions : a list containing the options for Alluvial plots

#' @details The list of parameters \code{plotOptions} provides the following options
#'  * "curvy" numeric, controls the curvature of the alluvial. Default value = 0.3
#'  * "alpha numeric, vector of transparency of the stripes. Default value = 0.8
#'  * "gap.width" numeric, relative width of inter-category gaps. Default value = 0.1
#'  * "col" vector of colors of the stripes. Default value is NULL: the stripes are colored according to the clusters of the first clustering
#'  * "border" vector of border colors for the stripes. Default is white
#'  * "reorder" boolean. If TRUE, the clusters of each clustering are stacked so as to follow the clusters of the previous one, which limits the crossings of the stripes. The labels of the clusters are not changed. Default value = TRUE
#'
#' @return display the alluvial plot, returns invisibly a list with the plotOptions and the table of frequencies
#' @export
#'
#' @examples
#' listMemberships <- list(C1 = rep(c('A','B','C'),each=10),C2 = rep(c(1,2,4),10))
#' plotAlluvial(listMemberships)

plotAlluvial = function(listMemberships,plotOptions = list()){

  nbMemb <- length(listMemberships)
  ### check sizes
  if (nbMemb < 2){stop('Provide at least two clusterings')}
  L <- sapply(listMemberships,length)
  test_size <- sum(abs(diff(L))) == 0
  if (!test_size){stop('The clusterings are not of equal sizes')}
  #-------------------------------------"
  currentOptions = list(
    curvy = 0.3,
    alpha = 0.8,
    gap.width = 0.1,
    col = NULL,
    border = "white",
    reorder = TRUE
  )
  currentOptions[names(plotOptions)] <- plotOptions
  #-------------------------------------"

  if (is.null(names(listMemberships))){names(listMemberships) = paste('C',c(1:nbMemb),sep='.')}
  U <- as.data.frame(lapply(listMemberships, as.factor), check.names = FALSE)

  # order the clusters of each clustering according to the mean position, in the
  # previous clustering, of their nodes (the first cluster is drawn at the bottom)
  if (currentOptions$reorder) {
    for (i in 2:nbMemb) {
      posPrev <- as.numeric(U[[i - 1]])
      score   <- tapply(posPrev, U[[i]], mean)
      U[[i]]  <- factor(U[[i]], levels = names(sort(score)))
    }
  }

  A <- as.data.frame(table(U))
  w   <- which(A$Freq!=0)
  A <- A[w,]
  if (nrow(A) == 1){stop('No clusterings to compare.')}

  # by default, one color per cluster of the first clustering
  col <- currentOptions$col
  if (is.null(col)) {
    palette <- grDevices::hcl.colors(nlevels(U[[1]]), "Dark 3")
    col <- palette[as.numeric(A[[1]])]
  }

  alluvial::alluvial(A[,1:nbMemb],freq=A$Freq , xw=currentOptions$curvy, alpha=currentOptions$alpha,
           gap.width=currentOptions$gap.width, col= col, border=currentOptions$border,
           axis_labels = names(listMemberships))
  invisible(list(plotOptions  = currentOptions, tableFreq  =A))
}
