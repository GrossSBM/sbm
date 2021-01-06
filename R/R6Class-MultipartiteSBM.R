#' R6 Class definition of a Multipartite SBM
#'
#' R6 virtual class for Multipartite SBM representation
#'
#' @import purrr dplyr
#' @export
MultipartiteSBM <-
  R6::R6Class(
    classname = "MultipartiteSBM",
    # fields for internal use (referring to the mathematical notation)
    private = list(
      model     = NULL, # list of characters describing the distributions of the edges (bernoulli, poisson, gaussian)
      listNet   = NULL, # list of SimpleSBMs and BipartiteSBMs composing the multipartite network
      E         = NULL, # matrix describing the organization of the multipartite network
      dimFG     = NULL, # number of nodes in each function groups
      namesFG   = NULL, # labels of the functional groups
      allZ      = NULL, # list of memberships
      pi        = NULL, # list of vectors of parameters for block prior probabilities
      theta     = NULL  # list of connectivity parameters between edges
    ),
    public = list(
      #' @description constructor for Multipartite SBM
      #' @param listSBM list of SimpleSBM or BipartiteSBM
      #' @param memberships list of memberships for each node in each functional group. Default value is NULL
      initialize = function(listSBM, memberships = NULL) {
        private$listNet <- listSBM
        private$namesFG <- listSBM %>% map("dimLabels") %>% unlist() %>% unique()

        # ###
        # E_FG <- lapply(listSBM,function(net){return(c(net$dimLabels$row,net$dimLabels$col))})
        # E_FG <- do.call(rbind,E_FG)
        # E <- matrix(sapply(E_FG,function(a){which(private$namesFG == a)}), self$nbNetworks,2)
        # private$E <-  E
        # private$dimFG <- sapply(1:self$nbLabels ,function(k){
        #   u <- which(E[,1] == k); v = 1;
        #   if (length(u) == 0) {u <- which(E[,2] == k); v = 2}
        #   u <- u[1]
        #   dim(listSBM[[u]]$netMatrix)[v]}
        # )

### alternative to above code with purrr
        private$dimFG <- listSBM %>% map("dimension") %>% unlist() %>% unique()
        private$E     <- listSBM %>% map_df("dimLabels") %>%
           map(factor, levels = private$namesFG) %>% map_df(as.numeric) %>% as.matrix()
###

        private$allZ <- memberships
        private$model <- map_chr(listSBM, "modelName")
      },
      #' @description print method
      #' @param type character to tune the displayed name
      show = function(type = "Multipartite Stochastic Block Model"){
        cat(type, "\n")
        cat(self$nbLabels, "functional groups (", self$dimLabels, "), ", self$nbNetworks, "networks\n")
        cat("=====================================================================\n")
        cat("nbNodes per FG = (", self$nbNodes, ") --  nbBlocks per FG = (",self$nbBlocks, ")\n")
        cat("distributions on each network =(", self$modelName ,")\n")
        cat("=====================================================================\n")
        cat("* Useful fields \n")
        cat(" $nbNetwork, $nbNodes, $nbBlocks, $dimLabels, $archiMultipartite \n")
        cat(" $modelName, $blockProp, $connectParam, $memberships, \n")
        cat("* Useful functions \n")
        cat("$plot, $optimize \n")
      },
      #' @description print method
      print = function() self$show(),
      #' @description plot Multipartite Network
      #' @param type character for the type of plot: either 'data' (true connection), 'expected' (fitted connection) or 'meso' (mesoscopic view). Default to 'data'.
      #' @param ordered TRUE is the matrices are plotted after reorganization with the blocks. Default value = TRUE
      #' @param plotOptions list of plot options for the mesoscopic view or matrix view
      plot = function(type = c('data','expected','meso'), ordered = TRUE, plotOptions = list()){

        type <- match.arg(type)

        if (type  == 'meso'){
          outP <-
            plotMesoMultipartite(
              private$E,
              private$theta,
              private$pi,
              private$model,
              self$directed,
              private$dimFG,
              private$namesFG,
              plotOptions
            )
        } else {
          listNetMatrix <- switch(type,
                                 'data' = map(private$listNet,"netMatrix"),
                                 'expected' = self$predict()
          )
          if (ordered) clust <- private$allZ else clust <- NULL
          outP <-
            plotMultipartiteMatrix(
              listNetMatrix,
              private$E,
              private$dimFG,
              private$namesFG,
              distrib  = private$model,
              clustering = clust,
              plotOptions = plotOptions
            )
        }
        outP
      }
    ),
    active = list(
      #' @field nbNetworks : number of networks in the multipartite network
      nbNetworks    = function(value) {length(private$listNet)},
      #' @field listSBM : list of SimpleSBMs or BipartiteSBMs
      listSBM    = function(value) {private$listNet},
      #' @field archiMultipartite : organization of the multipartite network
      archiMultipartite     = function(value) {private$E},
      #' @field dimLabels  : labels of the functional groups
      dimLabels   = function(value){private$namesFG},
      #' @field nbLabels  : number of Functional groups involved in the multipartite
      nbLabels   = function(value){length(private$namesFG)},
      #' @field nbNodes  : number of Nodes in each FG,
      nbNodes  = function(value){setNames(private$dimFG, private$namesFG)},
      #' @field expectation expected values of connection under the currently adjusted model
      expectation = function() {self$predict()},
      #' @field modelName vector of characters, the family of model for the distribution of the edges in each network
      modelName    = function(value) {private$model},
      #' @field directed : vector of boolean
      directed  = function(value){map(private$listNet, "directed") %>% map_lgl(~ifelse(is.null(.x), NA, .x))}
    )
  )

#' MultipartiteSBM Plot
#'
#' Basic matrix plot method for SBM object
#'
#' @param x an object inheriting from class MultipartiteSBM
#' @param type character for the type of plot: either 'data' (true connection) or 'expected' (fitted connection) or 'meso' (meso-scopic). Default to 'data'.
#' @param ordered logical: should the functional group be ordered according to the clustering? Default to \code{TRUE}.
#' @param plotOptions list with parameters.
#' @param ... additional parameters for S3 compatibility. Not used
#' @details The list of parameters \code{plotOptions} for the mesoscopic plot is
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
#' The list of parameters \code{plotOptions} for the matrix plot is
#' \itemize{
#'  \item{"normalized":}{Boolean. TRUE if the various matrices are presented in the same scale (between O and 1). FALSE otherwise. Default value FALSE}
#'  \item{"compact":}{Boolean. Default value is TRUE if you ask for the matrices to be transposed to have a more compact view}
#'  \item{"legend": }{Boolean. Set TRUE if you   want to see the legend. Default value is FALSE}
#'  \item{"legend.title": }{Boolean. Set TRUE if you want to print the title of the legend. Default value is FALSE}
#'  \item{"legend.position": }{Position of the legend. Possible values are 'bottom', 'top','left,'right'. Default value is 'bottom'}
#'  \item{"nodeNames": }{Set true if the node Names must be plotted. Default value is FALSE}
#'  \item{"line.color":}{The color of the lines to separate groups. Default value is red}
#'  \item{"line.width":}{Width  of the lines to separate groups. Default value is NULL, automatically chosen}
#'  \item{"title": }{Title of the plot. Default value is NULL}
#'  }
#' @return a ggplot2 object for the \code{'data'} and \code{'expected'}, a list with the igraph object \code{g}, the \code{layout} and the \code{plotOptions} for the \code{'meso'}
#' @export
plot.MultipartiteSBM = function(x, type = c('data', 'expected', 'meso'), ordered = TRUE, plotOptions = list(), ...){

  type <- match.arg(type)
  if (type == 'meso'){
    invisible(x$plot(type, ordered, plotOptions))
  } else {
    x$plot(type,ordered, plotOptions)
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

