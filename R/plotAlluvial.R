#' Plot an alluvial plot between  clusterings
#'
#' @param listMemberships  : a list vectors containing the memberships
#' @param plotOptions : a list containing the options for Alluvial plots

#' @details The list of parameters \code{plotOptions} provides the following options
#'  \itemize{
#'  \item{"curvy"}{numeric, controls the curvature of the alluvial. Default value = 0.3}
#'  \item{"alpha}{numeric, vector of transparency of the stripes. Default value = 0.8}
#'  \item{"gap.width"}{numeric, relative width of inter-category gaps. Default value = 0.1}
#'  \item{"col"}{vector of colors of the stripes. Default value = "darkolivegreen3"}
#'  \item{"border"}{vector of border colors for the stripes. Default is white}
#' }
#'
#' @return display the alluvial plot, returns the plotOptions as a list
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
    col = "darkolivegreen3",
    border="white"
  )
  currentOptions[names(plotOptions)] <- plotOptions
  #-------------------------------------"

  U <-as.data.frame(listMemberships[[1]])
  for (i in 1:nbMemb){
    U[,i] = listMemberships[[i]]
  }
  if (is.null(names(listMemberships))){names(listMemberships) = paste('C',c(1:nbMemb),sep='.')}
  names(U) <- names(listMemberships)
  A <- as.data.frame(table(U))
  w   <- which(A$Freq!=0)
  A <- A[w,]
  if (nrow(A) == 1){stop('No clusterings to compare.')}
  alluvial::alluvial(A[,1:nbMemb],freq=A$Freq , xw=currentOptions$curvy, alpha=currentOptions$alpha,
           gap.width=currentOptions$gap.width, col= currentOptions$col, border=currentOptions$border)
  return(list(plotOptions  = currentOptions, tableFreq  =A))
}
