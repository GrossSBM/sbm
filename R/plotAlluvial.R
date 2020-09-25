#' Plot an alluvial plot between  clusterings
#'
#' @param list  : a list vectors containing the memberships
#' @param plotOptions : a list containing the options for Alluvial plots
#' @return display the alluvial plot, return the plotOptions
#' @export
#'
#' @examples
#' listMemberships <- list(C1 = rep(c('A','B','C'),each=10),C2 = rep(c(1,2,4),10))
#' plotAlluvial(listMemberships)

plotAlluvial = function(listMemberships,plotOptions = list()){

  nbMemb <- length(listMemberships)
  ### check sizes
  if (nbMemb < 2){stop('Provide at least to clustergins')}
  L <- sapply(listMemberships,function(u){length(u)})
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
  alluvial::alluvial(A[,1:nbMemb],freq=A$Freq , xw=currentOptions$curvy, alpha=currentOptions$alpha,
           gap.width=currentOptions$gap.width, col= currentOptions$col, border=currentOptions$border)
  return(list(plotOptions  = currentOptions, tableFreq  =A))
}
