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
      #' @param blockProp parameters for block proportions (vector of list of vectors)
      #' @param connectParam list of parameters for connectivity
      #' @param covarParam optional vector of covariates effect
      #' @param covarList optional list of covariates data
      initialize = function(model='', dimension=numeric(2), blockProp=numeric(0), connectParam=list(mean = matrix()), covarParam=numeric(length(covarList)), covarList=list()) {

        ## SANITY CHECK
        stopifnot(is.character(model))
        stopifnot(model %in% available_models_edges)
        stopifnot(is.numeric(dimension), length(dimension) == 2)
        stopifnot(is.list(connectParam), is.matrix(connectParam$mean))
        stopifnot(all.equal(length(covarParam), length(covarList)))
        stopifnot(all(sapply(covarList, nrow) == dimension[1]))
        stopifnot(all(sapply(covarList, ncol) == dimension[2]))

        ## MODEL & PARAMETERS
        private$model <- model
        private$dim   <- dimension
        private$X     <- covarList
        private$pi    <- blockProp
        private$theta <- connectParam
        private$beta  <- covarParam
        private$link  <- switch(model,
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
      #' @param type character for the type of plot: either 'data' (true connection) or 'expected' (fitted connection). Default to 'data'.
      #' @param ordered logical: should the rows and columns be reordered according to the clustering? Default to \code{TRUE}.
      #' @param rowLabel character : type of the individual in row. Default to \code{NULL}.
      #' @param colLabel character : type of the individual in col. Default to \code{NULL}.
      #' @return a ggplot2 object
      #' @import ggplot2
      plot = function(type = c('data', 'expected'), ordered = TRUE, rowLabel = NULL, colLabel = NULL) {

        Mat <- switch(match.arg(type), data = self$netMatrix, expected = self$expectation)

        if (ordered) {
          if (is.vector(self$memberships)) {cl = list(row = self$memberships)}
          if (is.list(self$memberships)) {cl =  self$memberships; names(cl) = c('row','col') }
        }else{
          cl = NULL
        }

        if (is.null(colLabel)) {colLabel = ''}
        if (is.null(rowLabel)) {rowLabel = ''}

        P <- plotMatrix(Mat = Mat,rowFG = rowLabel,colFG = colLabel, clustering = cl)
        P
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
        cat("  $blockProp, $connectParam, covarParam, $covarList, $covarEffect \n")
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

## Auxiliary function to check the given class of an objet
is_SBM <- function(Robject) {inherits(Robject, "SBM")}

#' Extract model coefficients
#'
#' Extracts model coefficients from objects with class \code{\link[=SBM]{SBM}} and children (\code{\link[=SimpleSBM_fit]{SimpleSBM_fit}},
#' \code{\link[=BipartiteSBM_fit]{BipartiteSBM_fit}})
#'
#' @param object an R6 object inheriting from class SBM_fit (like SimpleSBM_fit or BipartiteSBM_fit)
#' @param type type of parameter that should be extracted. Either 'memberships' for \deqn{\pi}, 'connectivity' for \deqn{\theta},
#'  or "covariates" for \deqn{\beta}. Default is 'connectivity'.
#' @param ... additional parameters for S3 compatibility. Not used
#' @return vector or list of parameters.
#' @export
coef.SBM <- function(object, type = c( 'connectivity', 'membership', 'covariates'), ...) {
  stopifnot(is_SBM(object))
  switch(match.arg(type),
         membership   = object$blockProp,
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
#' @param rowLabel character : type of the individual in row. Default to \code{NULL}.
#' @param colLabel character : type of the individual in col. Default to \code{NULL}.
#' @param ... additional parameters for S3 compatibility. Not used
#' @return a ggplot2 object
#' @export
plot.SBM = function(x, type = c('data', 'expected'), ordered = TRUE, rowLabel = NULL, colLabel = NULL, ...){
  p <- x$plot(type, ordered, rowLabel, colLabel)
  p
}



