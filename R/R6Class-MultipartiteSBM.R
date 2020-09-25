#' R6 Class definition of a Multipartite SBM
#'
#' R6 virtual class for Multipartite SBM representation
#'
#' @export
MultipartiteSBM <-
  R6::R6Class(classname = "MultipartiteSBM",
     # fields for internal use (referring to the mathematical notation)
     private = list(
       nbNet = NULL,
       listNet = NULL,
       E = NULL,
       nbFG = NULL,
       dimFG = NULL,
       namesFG = NULL,
       allZ = NULL,
       directed  = NULL,
       distrib = NULL,
       pi  = NULL # list of vectors of parameters for block prior probabilities
       ),
     public = list(
       #' @description constructor for Multipartite SBM
       #' @param listSBM list of SimpleSBM or BipartiteSBM
       #' @param memberships list of memberships for each node in each function group.Default value is NULL
       initialize = function(listSBM, memberships = NULL) {
         private$listNet <- listSBM
         private$nbNet <- length(listSBM)
         private$namesFG <- unique(unlist(lapply(listSBM, function(net){net$dimLabels})))
         private$nbFG <- length(private$namesFG)
         E_FG <- lapply(listSBM,function(net){return(c(net$dimLabels$row,net$dimLabels$col))})
         E_FG <- do.call(rbind,E_FG)
         E <- matrix(sapply(E_FG,function(a){which(private$namesFG == a)}),private$nbNet,2,byrow = FALSE)
         private$E <-  E
         private$dimFG <- sapply(1:private$nbFG ,function(k){
           u <- which(E[,1] == k); v = 1;
           if (length(u) == 0) {u <- which(E[,2] == k); v = 2}
           u <- u[1]
           dim(listSBM[[u]]$netMatrix)[v]})
         private$allZ <- memberships
         private$distrib <- sapply(listSBM, function(net) {net$modelName})
         private$directed <- sapply(listSBM, function(net) {if(is.null(net$directed)){return(NA)}else{return(net$directed)}})
         },
         #' @description plot Multipartite Network
         #' @param type character for the type of plot: either 'data' (true connection), 'expected' (fitted connection) or 'meso' (mesoscopic view). Default to 'data'.
         #' @param normalized TRUE if the various matrices are renormalized. FALSE otherwise. Default value = FALSE
         #' @param ordered TRUE is the matrices are plotted after reorganization with the blocks. Default value = TRUE
         #' @param plotOptions list of plot options for the mesoscopic view
        plot = function(type=c('data','expected','meso'),normalized = FALSE, ordered = TRUE, plotOptions = list()){

          if (length(type)>1){type='data'}
          if (type %in% c('data','expected')){
            listNetMatrix = switch(type,
                                   'data'= lapply(private$listNet,function(s){s$netMatrix}),
                                   'expected' = self$predict()
                                   )
            if (ordered) { clust = private$allZ }else{ clust = NULL}
            g <- plotMultipartiteMatrix(listNetMatrix, private$E, private$dimFG, private$namesFG,normalized = normalized, clustering = clust)
            outP <- g
          }
          if (type  == 'meso'){
            layout <- NULL
            g <- NULL
            outP <- list(layout = layout,g=g )
          }
          outP
        }
     ),
     active = list(
         #' @field nbNetworks : number of networks in the multipartite network
         nbNetworks    = function(value) {private$nbNet},
         #' @field listSBM : list of SimpleSBMs or BipartiteSBMs
         listSBM    = function(value) {private$listNet},
         #' @field archiMultipartite : organization of the multipartite network
         archiMultipartite     = function(value) {private$E},
         #' @field dimLabels  : labels of the functional groups
         dimLabels   = function(value){private$namesFG},
         #' @field nbLabels  : number of Functional groups involved in the multipartite
         nbLabels   = function(value){private$nbFG},
         #' @field nbNodes  : number of Nodes in each FG,
         nbNodes  = function(value){private$dimFG},
         #' @field expectation expected values of connection under the currently adjusted model
         expectation = function() {self$predict()}

     )
     )



#' MultipartiteSBM Plot
#'
#' Basic matrix plot method for SBM object
#'
#' @param x an object inheriting from class MultipartiteSBM
#' @param type character for the type of plot: either 'data' (true connection) or 'expected' (fitted connection) or 'meso' (meso-scopic). Default to 'data'.
#' @param normalized TRUE if the various matrices are renormalized. FALSE otherwise. Default value = FALSE
#' @param ordered logical: should the functional group be ordered according to the clustering? Default to \code{TRUE}.
#' @param plotOptions list with parameters for 'meso' type plot
#' @param ... additional parameters for S3 compatibility. Not used
#' @details The list of parameters \code{plotOptions} is
#'  \itemize{
#'  \item{"seed": }{seed to control the layout}
#'  \item{"title": }{character string for the title. Default value is NULL}
#'  \item{"layout": }{Default value = NULL}
#'  \item{"vertex.color": }{Default value is "salmon2"}
#'  \item{"vertex.frame.color": }{Node border color.Default value is "black" }
#'  \item{"vertex.shape": }{One of "none", "circle", "square", "csquare", "rectangle" "crectangle", "vrectangle", "pie", "raster", or "sphere". Default value = "circle"}
#'  \item{"vertex.size": }{Size of the node (default is 2)}
#'  \item{"vertex.size2": }{The second size of the node (e.g. for a rectangle)}
#'  \item{"vertex.label": }{Names of the vertices. Default value is the label of the nodes}
#'  \item{"vertex.label.color": }{Default value is  "black"}
#'  \item{"vertex.label.font": }{Default value is 2. Font: 1 plain, 2 bold, 3, italic, 4 bold italic, 5 symbol}
#'  \item{"vertex.label.cex": }{Font size (multiplication factor, device-dependent).Default value is  0.9.}
#'  \item{"vertex.label.dist": }{Distance between the label and the vertex. Default value is  0}
#'  \item{"vertex.label.degree": }{The position of the label in relation to the vertex. default value is 0}
#'  \item{"edge.threshold": }{Threshold under which the edge is not plotted. Default value is = -Inf}
#'  \item{"edge.color": }{Default value is "gray"}
#'  \item{"edge.width": }{Factor parameter. Default value is 10}
#'  \item{"edge.arrow.size": }{Default value is 1}
#'  \item{"edge.arrow.width": }{Default value is 2}
#'  \item{"edge.lty": }{Line type, could be 0 or "blank", 1 or "solid", 2 or "dashed", 3 or "dotted", 4 or "dotdash", 5 or "longdash", 6 or "twodash". Default value is "solid"}
#'  \item{"edge.curved": }{Default value is = 0.3}
#' }
#' @return a ggplot2 object for the \code{'data'} and \code{'expected'}, a list with the igraph object \code{g}, the \code{layout} and the \code{plotOptions} for the \code{'meso'}
#' @export
plot.MultipartiteSBM = function(x, type = c('data', 'expected', 'meso'), normalized = FALSE, ordered = TRUE, plotOptions = list(), ...){

  if (length(type)>1){ type = 'data'}
  if (type=='meso'){
    invisible(x$plot(type, normalized, ordered, plotOptions))
  }else{
    x$plot(type,  normalized,ordered, plotOptions)
  }
}



#' Check  if an object is MultipartiteSBM
#'
#' Auxiliary function to check the given class of an object
#' @param  Robject an R6 object inheriting from class MultipartiteSBM
#' @return TRUE or FALSE
#' @export
is_MultipartiteSBM <- function(Robject) {inherits(Robject, " MultipartiteSBM")}


#' Model Predictions
#'
#' Make predictions from an Multipartite SBM.
#'
#' @param object an R6 object inheriting from class SBM_fit (like SimpleSBM_fit or BipartiteSBM_fit)
#' @param ... additional parameters for S3 compatibility. Not used
#' @return a list of matrices of expected values for each dyad
#' @export
predict.MultipartiteSBM <- function(object,...) {
  stopifnot(is_MultipartiteSBM(object))
  object$predict()
}

