#' Plot the matrices corresponding to a Multipartite Network
#'
#' @param listSBM  : a list of objects representing the multipartite network (see)
#' @param normalized TRUE if the various matrices are presented in the same scale (between O and 1). FALSE otherwise. Default value FALSE
#' @param memberships : a list of length equal to the number of Functional Groups providing clusterings inside these FGs
#' @param plotOptions : list().
#' @details plotOptions is a list containing the following items
#' \itemize{
#'  \item{"compact": }{Boolean. Default value is TRUE if you ask for the matrices to be transposed to have a more compact view}
#'  \item{"legend": }{Boolean. FALSE if you do not want to see the legend}
#'  \item{"line.color ": }{The color of the lines to separate groups. Default value is red}
#'  \item{"line.width ": }{Width  of the lines to separate groups. Default value is NULL, automatically chosen}
#'  }
#'
#' @return a ggplot object corresponding to the plot
#' @export
#'
#' @examples
#' data("multipartiteEcologicalNetwork")
#' Net <- multipartiteEcologicalNetwork
#' type='bipartite'
#' model = 'bernoulli'
#' directed = FALSE
#' PlantFlovis = defineSBM(Net$Inc_plant_flovis, model,type,directed,
#'                         dimLabels = list(row="Plants",col="Flovis"))
#' PlantAnt = defineSBM(Net$Inc_plant_ant,model,type,directed,
#'                      dimLabels =list(row = "Plants", col = "Ants"))
#' PlantBird = defineSBM(Net$Inc_plant_bird,model,type,directed,
#'                       dimLabels =list(row = "Plants",col = "Birds"))
#' plotMyMultipartiteMatrix(list(PlantFlovis,PlantAnt,PlantBird))
#'
#'

plotMyMultipartiteMatrix = function(listSBM,normalized = FALSE, memberships = NULL,plotOptions=list()){

  myMSBMObject <- MultipartiteSBM$new(listSBM,memberships = memberships)
  if(is.null(memberships)){ordered = FALSE}else{ordered = TRUE}
  g <- myMSBMObject$plot(type='data', normalized,ordered = ordered,plotOptions)
  g
}



