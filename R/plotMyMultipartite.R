#' Plot the matrices corresponding to a Multipartite Network
#'
#' @param listSBM  : a list of objects representing the multipartite network (see)
#' @param memberships : a list of length number of Functional Groups with clusterings inside these FG
#' @return a ggplot object corresponding to the plot
#' @export
#'
#' @examples
#' data(MPEcoNetwork, package = "GREMLIN")
#' PlantFlovis = defineSBM(MPEcoNetwork$Inc_plant_flovis, type = "bipartite",model = 'bernoulli',dimLabels = list(row="Plants",col="Flovis"))
#' PlantAnt = defineSBM(MPEcoNetwork$Inc_plant_ant,type = "bipartite",model = 'bernoulli',dimLabels =list(row = "Plants", col = "Ants"))
#' PlantBird = defineSBM(MPEcoNetwork$Inc_plant_bird,type = "bipartite",model = 'bernoulli',dimLabels =list(row = "Plants",col = "Birds"))
#' plotMyMultipartiteMatrix(list(PlantFlovis,PlantAnt,PlantBird))
#'
#'

plotMyMultipartiteMatrix = function(listSBM,normalizing = FALSE, memberships = NULL){

  myMSBMObject <- MultipartiteSBM$new(listSBM,memberships = memberships)
  g <- myMSBMObject$plot(normalizing)
  g
}

