#' R6 virtual class for SBM representation (mother class of Simple and Bipartite SBM fit and sampler)
#'
#' @import R6
#' @include R6Class-SBM.R
#' @export
SimpleSBM <- # this class inherit from SBM and allow to use multipartite as a list of SimpleSBM and BipartiteSBM
  R6::R6Class(classname = "SimpleSBM",
              inherit = SBM,
              private = list(
                directed_ =NULL,
                Y = NULL,
                tau = NULL

              ),
              public = list(
                #' @description constructor for a Simple SBM fit
                #' @param adjacencyMatrix square (weighted) matrix
                #' @param model character (\code{'bernoulli'}, \code{'poisson'}, \code{'gaussian'})
                #' @param directed logical, directed network or not. In not, \code{adjacencyMatrix} must be symmetric.
                #' @param dimLabels list of labels of each dimension (in row, in columns)
                #' @param covarList and optional list of covariates, each of whom must have the same dimension as \code{adjacencyMatrix}
                initialize = function(adjacencyMatrix, model, directed, dimLabels=list(row = "rowLabel", col = "colLabel"), covarList=list()) {

                  ## SANITY CHECKS
                  stopifnot(all.equal(nrow(adjacencyMatrix), ncol(adjacencyMatrix)))  # matrix must be square
                  stopifnot(isSymmetric(adjacencyMatrix) == !directed)                # symmetry and direction must agree

                  ## INITIALIZE THE SBM OBJECT ACCORDING TO THE DATA
                  super$initialize(model = model, dimension = dim(adjacencyMatrix), dimLabels = dimLabels, covarList = covarList)
                  private$Y <- adjacencyMatrix
                  private$directed_ <- directed
                },
                #' @description mesoscopic view for a Simple SBM
                #' @param plotOptions (see details in \code{plotMeso.SimpleSBM})
                plotMesoscopic = function(plotOptions){
                  if (is.null(private$theta)) {
                    stop('no mesoscopic plot possible. Estimate before')
                    }else{
                      plotOptions = list()
                      plotMeso(thetaMean = private$theta$mean,
                               pi = private$pi,
                               directed=private$directed_,
                               bipartite = FALSE,
                               nbNodes  = private$dimension[1],
                               nodeLabels = public$dimLabels$row,
                               plotOptions)
                      }
                }
              ),
              active = list(
                #' @field varProb variational probabilities of nodes being in a block
                varProb    = function(value) {if (missing(value)) return(private$tau) else {stopifnot(nrow(value)==private$dim[1])
                   private$tau <- value}},
                #' @field memberships vector of clustering
                memberships = function(value) {as_clustering(private$tau)},
                #' @field directed is the network directed or not
                directed = function(value) {private$directed_}
              )
              )

#' SimpleSBM  Mesoscopic Plot
#'
#' Mesoscopic view of SimpleSBM object
#' @description Mesoscopic view of SimpleSBM object
#' @param x a object inheriting from class simpleSBM
#' @param plotOptions a list a of options (see details)
#' @param ... additional parameters for S3 compatibility. Not used
#' @details The list of parameters \code{plotOptions} is
#'  \itemize{
#'  \item{"seed": }{seed to control the layout}
#'  \item{"title": }{chain of chararcter for the title. Default value is NULL}
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
#' @return
#' @export
plotMeso.SimpleSBM = function(x, plotOptions = list(), ...){
  x$plotMesoscopic(plotOptions)
}

