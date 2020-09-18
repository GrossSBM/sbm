available_models_edges <- c('bernoulli', 'poisson', 'gaussian')

#' R6 virtual class for SBM representation (mother class of Simple and Bipartite SBM fit and sampler)
#'
#' @import R6
SBM <- # this virtual class is the mother of all subtypes of SBM (Simple or Bipartite)
  R6::R6Class(classname = "SBM",
    ## fields for internal use (referring to the mathematical notation)
    private = list(
      model   = NULL, # characters, the model name: distribution of the edges (bernoulli, poisson, gaussian)
      link    = NULL, # the link function (GLM-like)
      invlink = NULL, # the inverse link function (GLM-like)
      dim     = NULL, # vector: number of nodes in row and in col
      dimlab  = NULL, # vector: the type of nodes in row and in col
      pi      = NULL, # vector of parameters for block prior probabilities
      theta   = NULL, # connectivity parameters between edges
      beta    = NULL, # vector of covariates parameters
      Y       = NULL, # data matrix  (dim[1] x dim[2])
      X       = NULL  # list of covariates (list of dim[1] x dim[2] matrices)
    ),
    public = list(
      #' @description constructor for SBM
      #' @param model character describing the type of model
      #' @param dimension dimension of the network matrix
      #' @param dimLabels labels of each dimension (in row, in columns)
      #' @param blockProp parameters for block proportions (vector of list of vectors)
      #' @param connectParam list of parameters for connectivity
      #' @param covarParam optional vector of covariates effect
      #' @param covarList optional list of covariates data
      initialize = function(model='', dimension=numeric(2), dimLabels=vector("list",2), blockProp=numeric(0), connectParam=list(mean = matrix()), covarParam=numeric(length(covarList)), covarList=list()) {

        ## SANITY CHECK
        stopifnot(is.character(model))
        stopifnot(model %in% available_models_edges)
        stopifnot(is.list(dimLabels), length(dimLabels) == 2)
        stopifnot(is.numeric(dimension), length(dimension) == 2)
        stopifnot(is.list(connectParam), is.matrix(connectParam$mean))
        stopifnot(all.equal(length(covarParam), length(covarList)))
        stopifnot(all(sapply(covarList, nrow) == dimension[1]))
        stopifnot(all(sapply(covarList, ncol) == dimension[2]))

        ## MODEL & PARAMETERS
        private$model   <- model
        private$dim     <- dimension
        private$dimlab  <- dimLabels
        names(private$dimlab) <- c("row","col") # names of dimlab
        private$X       <- covarList
        private$pi      <- blockProp
        private$theta   <- connectParam
        private$beta    <- covarParam
        private$link    <- switch(model,
                "gaussian"  = function(x) {x},
                "poisson"   = function(x) {log(x)},
                "bernoulli" = function(x) {.logit(x)},
                )
        private$invlink  <- switch(model,
                "gaussian"  = function(x) {x},
                "poisson"   = function(x) {exp(x)},
                "bernoulli" = function(x) {.logistic(x)},
                )
      },
      #' @description basic matrix plot method for SBM object
      #' @param type character for the type of plot: either 'data' (true connection), 'expected' (fitted connection) or 'meso' (mesoscopic view). Default to 'data'.
      #' @param ordered logical: should the rows and columns be reordered according to the clustering? Default to \code{TRUE}.
      #' @param plotOptions list with the parameters for meso plot (see details in \code{plotMeso.SimpleSBM}
      #' @return a ggplot2 object or a standard plot if mesoscopic view
      #' @import ggplot2
      plot = function(type = c('data','expected','meso'), ordered = TRUE,plotOptions = list()) {

        if (length(type) > 1){type='data'}
        if (type == 'meso'){
          if (is.vector(self$memberships)){bipartite = FALSE}
          if (is.list(self$memberships)) {bipartite = TRUE}
          plotMeso(thetaMean =  private$theta$mean,
                   pi = private$pi,
                   model = private$model,
                   directed= private$directed_,
                   bipartite = bipartite,
                   nbNodes  = self$dimension,
                   nodeLabels = self$dimLabels,
                   plotOptions)
        }else{
            Mat <- switch(match.arg(type), data = self$netMatrix, expected = self$expectation)
            if (ordered) {
              if (is.vector(self$memberships)) {
                cl = list(row = self$memberships)
              }
              if (is.list(self$memberships)) {
                cl =  self$memberships
                names(cl) = c('row', 'col')
              }
            } else {
              cl = NULL
            }
          P <- plotMatrix(Mat = Mat, dimLabels = self$dimLabels, clustering = cl)
          P
          }

        },
      #' @description print method
      #' @param type character to tune the displayed name
      show = function(type = "Stochastic Block Model") {
        cat(type, "--", self$modelName, "variant\n")
        cat("=====================================================================\n")
        cat("Dimension = (", self$dimension, ") - (",
            self$nbBlocks, ") blocks and",
          ifelse(self$nbCovariates > 0, self$nbCovariates, "no"), "covariate(s).\n")
        cat("=====================================================================\n")
        cat("* Useful fields \n")
        cat("  $dimension, $modelName, $nbNodes, $nbBlocks, $nbCovariates, $nbDyads\n")
        cat("  $blockProp, $connectParam, covarParam, $covarList, $covarEffect, $dimLabels \n")
      },
      #' @description print method
      print = function() self$show()
    ),
    ## active binding to access fields outside the class
    active = list(
      #' @field dimension size-2 vector: dimension of the network
      dimension    = function(value) {private$dim},
      #' @field modelName character, the family of model for the distribution of the edges
      modelName    = function(value) {private$model},
      #' @field dimLabels vector of characters, the label of each dimension
      dimLabels    = function(value) {private$dimlab},
      #' @field nbCovariates integer, the number of covariates
      nbCovariates = function(value) {length(private$X)},
      #' @field blockProp vector of block proportions (aka prior probabilities of each block)
      blockProp    = function(value) {if (missing(value)) return(private$pi) else private$pi <- value},
      #' @field connectParam parameters associated to the connectivity of the SBM, e.g. matrix of inter/inter block probabilities when model is Bernoulli
      connectParam = function(value) {if (missing(value)) return(private$theta) else private$theta <- value},
      #' @field covarParam vector of regression parameters associated with the covariates.
      covarParam   = function(value) {if (missing(value)) return(private$beta) else private$beta <- value},
      #' @field covarList list of matrices of covariates
      covarList    = function(value) {if (missing(value)) return(private$X) else private$X <- value},
      #' @field covarEffect effect of covariates
      covarEffect  = function(value) {if (self$nbCovariates > 0) return(roundProduct(simplify2array(private$X), private$beta)) else return(numeric(0))},
      #' @field netMatrix the matrix (adjacency or incidence) encoding the network
      netMatrix    = function(value) {if (missing(value)) return(private$Y) else private$Y <- value}
    )
  )


# ========================================================================================
# PUBLIC S3 METHODS FOR SBM

#' Auxiliary function to check the given class of an object
#'
#' Auxiliary function to check the given class of an objet
#' @param  object an R6 object inheriting from class SBM
#' @return TRUE or FALSE
#' @export
is_SBM <- function(Robject) {inherits(Robject, "SBM")}

#' Extract model coefficients
#'
#' Extracts model coefficients from objects with class \code{\link[=SBM]{SBM}} and children (\code{\link[=SimpleSBM_fit]{SimpleSBM_fit}},
#' \code{\link[=BipartiteSBM_fit]{BipartiteSBM_fit}})
#'
#' @param object an R6 object inheriting from class SBM_fit (like SimpleSBM_fit or BipartiteSBM_fit)
#' @param type type of parameter that should be extracted. Either 'block' for \deqn{\pi}, 'connectivity' for \deqn{\theta},
#'  or "covariates" for \deqn{\beta}. Default is 'connectivity'.
#' @param ... additional parameters for S3 compatibility. Not used
#' @return vector or list of parameters.
#' @export
coef.SBM <- function(object, type = c( 'connectivity', 'block', 'covariates'), ...) {
  stopifnot(is_SBM(object))
  switch(match.arg(type),
         block        = object$blockProp,
         connectivity = object$connectParam,
         covariates   = object$covarParam)
}

#' Model Predictions
#'
#' Make predictions from an SBM.
#'
#' @param object an R6 object inheriting from class SBM_fit (like SimpleSBM_fit or BipartiteSBM_fit)
#' @param covarList a list of covariates. By default, we use the covariates associated with the model.
#' @param ... additional parameters for S3 compatibility. Not used
#' @return a matrix of expected values for each dyad
#' @importFrom stats predict
#' @export
predict.SBM <- function(object, covarList = object$covarList, ...) {
  stopifnot(is_SBM(object))
  object$predict(covarList)
}

#' SBM Plot
#'
#' Basic matrix plot method for SBM object
#' @description basic matrix plot method for SBM object
#' @param x a object inheriting from class SBM
#' @param type character for the type of plot: either 'data' (true connection) or 'expected' (fitted connection). Default to 'data'.
#' @param ordered logical: should the rows and columns be reoredered according to the clustering? Default to \code{TRUE}.
#' @param plotOptions list  with parameters for 'meso' type plot
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
#' @return a ggplot2 object of a standard plot for 'meso' plot
#' @export
plot.SBM = function(x, type = c('data', 'expected','meso'), ordered = TRUE, plotOptions=list(), ...){

  x$plot(type, ordered,plotOptions)

}



