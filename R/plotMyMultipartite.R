#' Plot the matrices corresponding to a Multipartite Network
#'
#' @param listSBM  : a list of objects representing the multipartite network (see)
#' @param normalized TRUE if the various matrices are presented in the same scale (between O and 1). FALSE otherwise. Default value FALSE
#' @param memberships : a list of length number of Functional Groups with clusterings inside these FG
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

plotMyMultipartiteMatrix = function(listSBM,normalized = FALSE, memberships = NULL){

  myMSBMObject <- MultipartiteSBM$new(listSBM,memberships = memberships)
  g <- myMSBMObject$plot(type='data', normalized)
  g
}



