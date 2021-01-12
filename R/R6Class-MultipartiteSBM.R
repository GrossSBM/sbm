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
      model     = NULL, # vector of characters for the model name: distributions of the edges
      directed_ = NULL, # vector of logical indicating if networks are directed, when appropriate
      netList   = NULL, # list of SimpleSBMs and BipartiteSBMs composing the multipartite network
      arch      = NULL, # matrix describing the organization of the multipartite network
      dim       = NULL, # number of nodes in each function groups
      dimlab    = NULL, # labels of the functional groups
      pi        = NULL, # list of vectors of parameters for block prior probabilities
      theta     = NULL  # list of connectivity parameters between edges
    ),
    public = list(
      #' @description constructor for Multipartite SBM
      #' @param model character describing the type of model
      #' @param architecture a 2-column matrix describing interactions between the networks
      #' @param directed vector of logical: are the network directed or not?
      #' @param dimension number of nodes in each functional groups
      #' @param dimLabels labels of each functional groups
      #' @param blockProp parameters for block proportions (vector of list of vectors)
      #' @param connectParam parameters of connectivity (vector of list of vectors)
      initialize = function(model = character(0), architecture = matrix(NA, 0, 2), directed = logical(0),
                            dimension = numeric(0), dimLabels = character(0), blockProp=list(), connectParam=list()) {

        ## SANITY CHECK
        stopifnot(is.character(model), model %in% available_models_edges)
        stopifnot(is.matrix(architecture), ncol(architecture) == 2)
        stopifnot(is.logical(directed), nrow(architecture) == length(directed))
        stopifnot(is.character(dimLabels), length(dimLabels) == length(dimension))
        stopifnot(is.list(connectParam))

        ## MODEL & PARAMETERS
        private$model  <- model
        private$arch   <- architecture
        private$dim    <- dimension
        private$dimlab <- dimLabels
        private$pi     <- blockProp
        private$theta  <- connectParam
        private$directed_ <- directed

      },
      #' @description print method
      #' @param type character to tune the displayed name
      show = function(type = "Multipartite Stochastic Block Model"){
        cat(type, "\n")
        cat(self$nbLabels, "functional groups (", self$dimLabels, "), ", self$nbNetworks, "networks\n")
        cat("=====================================================================\n")
        cat("nbNodes per FG = (", self$nbNodes, ") --  nbBlocks per FG = (",self$nbBlocks, ")\n")
        cat("distributions on each network: ", self$modelName ,"\n")
        cat("=====================================================================\n")
        cat("* Useful fields \n")
        cat("  $nbNetwork, $nbNodes, $nbBlocks, $dimLabels, $architecture \n")
        cat("  $modelName, $blockProp, $connectParam, $memberships, $networkList\n")
      },
      #' @description print method
      print = function() self$show(),
      #' @description plot Multipartite Network
      #' @param type character for the type of plot: either 'data' (true connection), 'expected' (fitted connection) or 'meso' (mesoscopic view). Default to 'data'.
      #' @param ordered TRUE is the matrices are plotted after reorganization with the blocks. Default value = TRUE
      #' @param plotOptions list of plot options for the mesoscopic view or matrix view
      plot = function(type = c('data','expected','meso'), ordered = TRUE, plotOptions = list()){

        if (ordered) clustering <- self$memberships else clustering <- NULL

        switch(match.arg(type),
          "meso" =
            plotMesoMultipartite(
              private$arch, self$connectParam, private$pi, private$model,
              private$directed_, private$dim, private$dimlab, plotOptions
            ),
          "data" =
            plotMultipartiteMatrix(
              map(private$netList,"netMatrix"),
              private$arch, private$dim, private$dimlab,
              private$model, clustering, plotOptions
            ),
          "expected" =
            plotMultipartiteMatrix(
              self$predict(),
              private$arch, private$dim, private$dimlab,
              private$model, clustering, plotOptions
            )
        )
      }
    ),
    active = list(
      #' @field modelName vector of characters, the family of model for the distribution of the edges in each network
      modelName    = function(value) {private$model},
      #' @field architecture organization of the multipartite network
      architecture = function(value) {private$arch},
      #' @field nbNetworks number of networks in the multipartite network
      nbNetworks = function(value) {length(private$directed_)},
      #' @field directed vector of boolean
      directed = function(value){private$directed_},
      #' @field dimension number of Nodes in each functional group,
      dimension = function(value){setNames(private$dim, private$dimlab)},
      #' @field nbNodes number of Nodes in each functional group
      nbNodes = function(value){setNames(private$dim, private$dimlab)},
      #' @field dimLabels labels of the functional groups
      dimLabels = function(value){private$dimlab},
      #' @field nbLabels number of functional groups involved in the multipartite
      nbLabels  = function(value){length(private$dimlab)},
      #' @field blockProp  block proportions in each function group
      blockProp = function(value) {private$pi},
      #' @field connectParam connection parameters in each network
      connectParam = function(value) {private$theta},
      #' @field networkList list of SimpleSBMs or BipartiteSBMs
      networkList = function(value) {private$netList},
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
#'  \item{"vertex.size": }{Size of the nodes (default factor is 1). Vector of length the number of FG}
#'  \item{"vertex.size2": }{The second size of the node (e.g. for a rectangle)}
#'  \item{"vertex.label": }{Names of the vertices. Default value is the label of the nodes}
#'  \item{"vertex.label.color": }{Default value is  "black"}
#'  \item{"vertex.label.font": }{Default value is 2. Font: 1 plain, 2 bold, 3, italic, 4 bold italic, 5 symbol}
#'  \item{"vertex.label.cex": }{Font size (multiplication factor, device-dependent).Default value is  0.9.}
#'  \item{"vertex.label.dist": }{Distance between the label and the vertex. Default value is  0}
#'  \item{"vertex.label.degree": }{The position of the label in relation to the vertex. default value is 0}
#'  \item{"edge.threshold": }{Threshold under which the edge is not plotted. Default value is = -Inf}
#'  \item{"edge.color": }{Default value is "gray"}
#'  \item{"edge.width": }{Factor parameter. Default value is 1}
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
    x$plot(type, ordered, plotOptions)
  }
}

#' Check  if an object is MultipartiteSBM
#'
#' Auxiliary function to check the given class of an object
#' @param  Robject an R6 object inheriting from class MultipartiteSBM
#' @return TRUE or FALSE
#' @export
is_MultipartiteSBM <- function(Robject) {inherits(Robject,"MultipartiteSBM")}

#' Model Predictions
#'
#' Make predictions from an Multipartite SBM.
#'
#' @param object an R6 object inheriting from class SBM_fit (like SimpleSBM_fit or BipartiteSBM_fit)
#' @param ... additional parameters for S3 compatibility. Not used
#' @return a list of matrices of expected values for each dyad
#' @export
predict.MultipartiteSBM <- function(object, ...) {
  stopifnot(is_MultipartiteSBM(object))
  object$predict()
}

#' Extract model coefficients
#'
#' Extracts model coefficients from objects with class \code{\link[=MultipartiteSBM]{MultipartiteSBM}} and children
#'
#' @param object an R6 object inheriting from class MultipartiteSBM
#' @param type type of parameter that should be extracted. Either 'block' for \deqn{\pi}, 'connectivity' for \deqn{\theta},
#'  or "covariates" for \deqn{\beta}. Default is 'connectivity'.
#' @param ... additional parameters for S3 compatibility. Not used
#' @return vector or list of parameters.
#' @export
coef.MultipartiteSBM <- function(object, type = c( 'connectivity', 'block'), ...) {
  stopifnot(is_MultipartiteSBM(object))
  switch(match.arg(type),
         block        = object$blockProp,
         connectivity = object$connectParam)
}
