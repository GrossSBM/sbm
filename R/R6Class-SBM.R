available_models_edges <- c('bernoulli', 'poisson', 'gaussian','ZIgaussian')

#' R6 virtual class for SBM representation (mother class of SimpleSBM, BipartiteSBM, MultipartiteSBM)
#'
#' @import R6
SBM <- # this virtual class is the mother of all subtypes of SBM (Simple or Bipartite)
  R6::R6Class(classname = "SBM",
    ## fields for internal use (referring to the mathematical notation)
    private = list(
      model         = NULL, # characters, the model name: distribution of the edges (bernoulli, poisson, gaussian)
      directed_     = NULL, # vector of logical indicating if networks are directed, when appropriate
      link          = NULL, # the link function (GLM-like)
      invlink       = NULL, # the inverse link function (GLM-like)
      dim           = NULL, # dimension: number of nodes for each group
      dimlab        = NULL, # vector: the type of nodes in row and in col
      pi            = NULL, # vector of parameters for block prior probabilities
      theta         = NULL, # connectivity parameters between edges
      beta          = NULL, # vector of covariates parameters
      Y             = NULL, # network data (matrix or list of matrices)
      X             = NULL, # list of covariates
      Z             = NULL, # indicator/probablities of blocks belonging
      sampling_func = NULL  # a list of functions to sample edge values, depending on the model
    ),
    public = list(
      #' @description constructor for SBM
      #' @param model character describing the type of model
      #' @param directed logical describing if the network data is directed or not
      #' @param dimension dimension of the network data
      #' @param dimLabels labels of each dimension
      #' @param blockProp parameters for block proportions (vector or list of vectors)
      #' @param connectParam list of parameters for connectivity
      #' @param covarParam optional vector of covariates effect
      #' @param covarList optional list of covariates data
      initialize = function(model        = vector("character", 0),
                            directed     = vector("logical"  , 0),
                            dimension    = vector("numeric"  , 0),
                            dimLabels    = vector("character", 0),
                            blockProp    = vector("numeric"  , 0),
                            connectParam = vector("list"     , 0),
                            covarParam   = numeric(length(covarList)),
                            covarList    = list()) {
        ## SANITY CHECK
        stopifnot(is.character(model), all(model %in% available_models_edges))
        stopifnot(is.logical(directed))
        stopifnot(is.character(dimLabels), length(dimLabels) == length(dimension))
        stopifnot(is.numeric(dimension), all(dimension > 0))
        stopifnot(is.list(connectParam))
        stopifnot(all.equal(length(covarParam), length(covarList)))

        ## MODEL & PARAMETERS
        private$model      <- model
        private$directed_  <- directed
        private$dim        <- dimension
        private$dimlab     <- dimLabels
        private$X          <- covarList
        private$pi         <- blockProp
        private$theta      <- connectParam
        private$beta       <- covarParam

        private$link <- map(model,
          ~switch(.x,
                "gaussian"   = function(x) {x},
                "ZIgaussian" = function(x) {x},
                "poisson"    = function(x) {log(x)},
                "bernoulli"  = function(x) {.logit(x)},
                )
          )
        private$invlink <- map(model,
          ~switch(.x,
                "gaussian"   = function(x) {x},
                "ZIgaussian" = function(x) {x},
                "poisson"    = function(x) {exp(x)},
                "bernoulli"  = function(x) {.logistic(x)},
                )
          )

        private$sampling_func <- map(model,
          ~switch(.x,
            "gaussian"   = function(n, param) rnorm (n = n, mean = param$mean, sd = sqrt(param$var)),
            "ZIgaussian" = function(n, param) rbinom(n = n, size = 1, prob = 1 - param$p0) * rnorm(n = n, param$mean, sd = sqrt(param$var)),
            "poisson"    = function(n, param) rpois (n = n, lambda = param$mean) ,
            "bernoulli"  = function(n, param) rbinom(n = n, size = 1, prob   = param$mean)
          )
        )
      },
      #' @description a method to sample a network data for the current SBM (blocks and edges)
      #' @param store should the sampled network be stored (and overwrite the existing data)? Default to FALSE
      #' @return a list with the sampled block and network
      rNetwork = function(store = FALSE) {
        Z <- self$rMemberships(store = store)
        E <- self$rEdges(store = store)
        list(indMemberships = Z, networkData = E)
      },
      #' @description print method
      #' @param type character to tune the displayed name
      show = function(type = "Stochastic Block Model") {
        cat(type, "--", self$modelName, "variant\n")
        cat("=====================================================================\n")
        cat("Dimension = (", self$nbNodes, ") - (",
            self$nbBlocks, ") blocks and",
          ifelse(self$nbCovariates > 0, self$nbCovariates, "no"), "covariate(s).\n")
        cat("=====================================================================\n")
        cat("* Useful fields \n")
        cat("  $nbNodes, $modelName, $dimLabels, $nbBlocks, $nbCovariates, $nbDyads\n")
        cat("  $blockProp, $connectParam, $covarParam, $covarList, $covarEffect \n")
        cat("  $expectation, $indMemberships, $memberships \n")
        cat("* R6 and S3 methods \n")
        cat("  $rNetwork, $rMemberships, $rEdges, plot, print, coef \n")
        },
      #' @description print method
      print = function() self$show()
    ),
    ## active binding to access fields outside the class
    active = list(
      #' @field modelName character, the family of model for the distribution of the edges
      modelName    = function(value) {private$model},
      #' @field directed mode of the network data (directed or not or not applicable)
      directed = function(value) {private$directed_},
      #' @field dimLabels vector or list of characters, the label of each dimension
      dimLabels    = function(value) {private$dimlab},
      #' @field nbNodes vector describing the number of the successive elements connecting the network
      nbNodes = function(value) {setNames(private$dim, private$dimlab)},
      #' @field nbCovariates integer, the number of covariates
      nbCovariates = function(value) {length(private$X)},
      #' @field blockProp block proportions (aka prior probabilities of each block)
      blockProp   = function(value) {private$pi},
      #' @field connectParam parameters associated to the connectivity of the SBM, e.g. matrix of inter/inter block probabilities when model is Bernoulli
      connectParam = function(value) {private$theta},
      #' @field covarParam vector of regression parameters associated with the covariates.
      covarParam  = function(value) {private$beta},
      #' @field covarList list of matrices of covariates
      covarList   = function(value) {private$X},
      #' @field covarArray the array of covariates
      covarArray  = function(value) {if (self$nbCovariates > 0) simplify2array(private$X) else return(array())},
      #' @field covarEffect effect of covariates
      covarEffect = function(value) {if (self$nbCovariates > 0) return(roundProduct(private$X, private$beta)) else return(numeric(0))},
      #' @field networkData the network data (adjacency or incidence matrix or list of such object)
      networkData = function(value) {return(private$Y)},
      #' @field expectation expected values of connection under the current model
      expectation = function() {self$predict()}
    )
  )


# ========================================================================================
# PUBLIC S3 METHODS FOR SBM

#' Auxiliary function to check the given class of an object
#'
#' Auxiliary function to check the given class of an object
#' @param  Robject an R6 object inheriting from class SBM
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
#' @param theta_p0 a threshold...
#' @param ... additional parameters for S3 compatibility. Not used
#' @return a matrix of expected values for each dyad
#' @importFrom stats predict
#' @export
predict.SBM <- function(object, covarList = object$covarList, theta_p0 = 0, ...) {
  stopifnot(is_SBM(object))
  object$predict(covarList, theta_p0)
}

#' SBM Plot
#'
#' Basic matrix plot method for SBM object or mesoscopic view
#'
#' @param x an object inheriting from class SBM
#' @param type character for the type of plot: either 'data' (true connection),  'expected' (fitted connection) or 'meso' (mesoscopic). Default to 'data'.
#' @param ordered logical: should the rows and columns be ordered according to the clustering? Default to \code{TRUE} (not taken into account for 'meso').
#' @param plotOptions list with parameters for 'meso' type plot and data type plot. Details are given below
#' @param ... additional parameters for S3 compatibility. Not used
#' @details The list of parameters \code{plotOptions}  for the mesoscopic plot is:
#'  \itemize{
#'  \item{"seed": }{seed to control the layout}
#'  \item{"title": }{character string for the title. Default value is NULL}
#'  \item{"layout": }{Default value = NULL}
#'  \item{"vertex.color": }{Default value is "salmon2"}
#'  \item{"vertex.frame.color": }{Node border color.Default value is "black" }
#'  \item{"vertex.shape": }{One of "none", "circle", "square", "csquare", "rectangle" "crectangle", "vrectangle", "pie", "raster", or "sphere". Default value = "circle"}
#'  \item{"vertex.size": }{Size of the node (default is 2)}
#'  \item{"vertex.size2": }{The second size of the node (e.g. for a rectangle)}
#'  \item{"vertex.label.name": }{Names of the vertices. Default value is the label of the nodes}
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
#'  \item{"edge.curved": }{Default value is = 0.3.}
#' }
#' For type = 'data' or 'expected plot', the list of parameters \code{plotOptions} is
#' \itemize{
#'  \item{"legend": }{Boolean. Set TRUE if you want to see the legend. Default value is FALSE}
#'  \item{"legend.title":}{Boolean. Set TRUE if you want to print the title of the legend. Default value is FALSE}
#'  \item{"legend.position":}{Position of the legend. Possible values are 'bottom', 'top','left,'right'. Default value is 'bottom'}
#'  \item{"rowNames":}{Set true if the rownames must be plotted. Default value is FALSE}
#'  \item{"colNames":}{Set true if the colNames must be plotted. Default value is FALSE}
#'  \item{"line.color": }{Chain of character. The color of the lines to separate groups if a clustering is provided. Default value is red}
#'  \item{"line.width": }{Numeric. Width  of the lines to separate groups. Default value is NULL, automatically chosen}
#'  \item{"title": }{Chain of character. Title of the plot. Default value is NULL}
#'  }
#' @return a ggplot2 object for the \code{'data'} and \code{'expected'}, a list with the igraph object \code{g} and the \code{layout} for the \code{'meso'}
#' @export
plot.SBM = function(x, type = c('data', 'expected', 'meso'), ordered = TRUE, plotOptions = list(), ...){

  stopifnot(is_SBM(x))
  type <- match.arg(type)
  if (type == 'meso'){
    invisible(x$plot(type, ordered, plotOptions))
  } else {
    x$plot(type, ordered, plotOptions)
  }
}

#' Extract model fitted values
#'
#' Extracts fitted values for object with class (\code{\link[=SimpleSBM_fit]{SimpleSBM_fit}},
#' \code{\link[=BipartiteSBM_fit]{BipartiteSBM_fit}}) or \code{\link[=MultipartiteSBM_fit]{multipartitepartiteSBM_fit}})
#'
#' @param object an R6 object inheriting from SimpleSBM_fit,  BipartiteSBM_fit or MultipartiteSBM_fit
#' @param ... additional parameters for S3 compatibility. Not used
#' @return a matrix of expected fitted values for each dyad
#' @importFrom stats fitted
#' @export
fitted.SBM <- function(object,  ...) {
  stopifnot(is_SBM(object))
  stopifnot(inherits(object, "SimpleSBM") |
            inherits(object, "BipartiteSBM") |
            inherits(object, "MultipartiteSBM")
      )
  object$predict()
}
