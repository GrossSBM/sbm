
#' R6 virtual class for ScoreSBM representation
#'
#' @import R6
ScoreSBM <- # this virtual class is the mother of all subtypes of SBM (Simple or Bipartite)
  R6::R6Class(classname = "ScoreSBM",
              ## fields for internal use (referring to the mathematical notation)
              private = list(
                model    = NULL, # distribution of the scores (gaussian by default )
                dim      = NULL, # nbNodes of the underlying network
                d        = NULL, # nb de Scores
                directed = NULL, # TRUE or FALSE
                theta    = NULL, # parameters of the model
                Y        = NULL  # list of Score Matrices data matrix  (dim[1] x dim[2])
              ),
              public = list(
                #' @description constructor for SBM
                #' @param model character describing the type of model
                #' @param nbScores number of scores
                #' @param dimension dimension of the network matrix (nb of nodes)
                #' @param directed directed (true or false)
                #' @param blockProp parameters for block proportions (vector of list of vectors)
                #' @param connectParam matrix of parameters for connectivity
                #' @param emissionParam parameters of the distribution of the Scores
                initialize = function(model = NA, nbScores = NA, dimension = NA, directed = NA, blockProp=NA, connectParam = NA, emissionParam = NA) {

                  ## MODEL & PARAMETERS
                  private$model <- model
                  private$d <- nbScores
                  private$dim   <- dimension
                  private$theta <- list(blockProp = blockProp)
                  private$theta$connectParam <- connectParam
                  private$theta$emissionParam <- emissionParam
                },
                #' @description basic matrix plot method for Noisy SBM object
                #' @param type character for the type of plot: either 'data' (true connection) or 'predictedNetwork' (fitted connection). Default to 'data'.
                #' @param ordered logical: should the rows and columns be reoredered accordigin to the clustering? Default to \code{TRUE}.
                #' @param color logical. Adapt colormap to the clustering. Default to \code{TRUE}.
                #' @importFrom corrplot corrplot
                plot = function(type = c('data', 'predictedNetwork'), ordered = TRUE, color = TRUE) {
                  mat <- switch(match.arg(type), data = self$listScore, predictedNetwork = list(self$predictedG()))
                  d <- length(mat)
                  cl <- self$memberships
                  Z <- as_indicator(as.factor(cl))
                  colors <- matrix(-ncol(Z), ncol(Z), ncol(Z)); diag(colors) <- floor(ncol(Z)/2) + (1:ncol(Z)) # discriminate intra/inter cols
                  colorMat <- Z %*% colors %*% t(Z)
                  if (ordered) {
                    ocl <- order(cl)
                    colorMat <- colorMat[ocl,ocl]
                    mat <- lapply(1:d,function(i){mat[[i]][ocl, ocl]})
                  }

                  if  (d <= 3) {sizeplot <- c(1,d)} else {sizeplot = c(floor(d/3) + 1,3)}
                  par(mfrow = sizeplot)
                  for (i in 1:d) {
                      corrplot(mat[[d]] * colorMat, is.corr = FALSE, tl.pos = "n", method = "color", cl.pos = "n", mar = c(0,0,1,0))
                  }
                },
                #' @description print method
                #' @param type character to tune the displayed name
                show = function(type = "Score Stochastic Block Model") {
                  cat(type, "--", self$modelName, "variant\n")
                  cat("=====================================================================\n")
                  cat("Dimension = (", self$dimension, ") - (",
                      self$nbBlocks, ") blocks and", self$nbScores, "Score(s)).\n")
                  cat("=====================================================================\n")
                  cat("* Useful fields \n")
                  cat("  $dimension, $modelName, $nbNodes, $nbBlocks, $nbScores, $nbDyads\n")
                  cat("  $blockProp, $connectParam, $emissionParam \n")
                },
                #' @description print method
                print =  function() self$show()
              ),
              ## active binding to access fields outside the class
              active = list(
                #' @field dimension size-2 vector: dimension of the network
                dimension       = function(value) {private$dim},
                #' @field modelName character, the family of model for the distribution of the edges
                modelName    = function(value) {private$model},
                #' @field blockProp vector of block proportions (aka prior probabilities of each block)
                blockProp     = function(value) {if (missing(value)) return(private$theta$blockProp)     else private$theta$blockProp     <- value},
                #' @field connectParam parameters associated to the connectivity of the SBM, e.g. matrix of inter/inter block probabilities when model is Bernoulli
                connectParam = function(value) {if (missing(value)) return(private$theta$connectParam)  else private$theta$connectParam  <- value},
                #' @field emissionParam emission parameters for the Scores.
                emissionParam   = function(value) {if (missing(value)) return(private$theta$emissionParam)   else private$theta$emissionParam  <- value},
                #' @field listScore the matrix (adjacency or incidence) encoding the network
                listScore    = function(value) {if (missing(value)) return(private$Y)      else private$Y      <- listScore}
              )
  )



# ScoreSBM$set("public", "sample", function(x = 1){
#
# }
#
# )




# ========================================================================================
# PUBLIC S3 METHODS FOR SBM

## Auxiliary function to check the given class of an objet
is_ScoreSBM <- function(Robject) {inherits(Robject, "ScoreSBM")}

#' Extract model coefficients
#'
#' Extracts model coefficients from objects with class \code{\link[=ScoreSBM]{ScoreSBM}}
#'
#' @param object an R6 object inheriting from class ScoreSBM_fit
#' @param type type of parameter that should be extracted. Either 'memberships' for \deqn{\pi}, 'connectivity' for \deqn{\alpha},
#'  or "emission" for \deqn{\mu_0,\Sigma_0, \mu_1, \sigma_1}.
#' @param ... additional parameters for S3 compatibility. Not used
#' @return vector or list of parameters.
#' @export
coef.ScoreSBM <- function(object, type = c( 'connectivity', 'membership', 'emission'), ...) {
  stopifnot(is_SBM(object))
  switch(match.arg(type),
         membership   = object$blockProp,
         connectivity = object$connectParam,
         emission   = object$emissionParam)
}

#' ScoreSBM Plot
#'
#' Basic matrix plot method for ScoreSBM object
#' @param x an R6 object inheriting from class SBM_fit (like SimpleSBM_fit or BipartiteSBM_fit)
#' @param type character for the type of plot: either 'data' (true connection) or 'estimatedNetwork' (fitted connection). Default to 'data'.
#' @param ordered logical: should the rows and columns be reordered according to the clustering? Default to \code{TRUE}.
#' @param color logical. Adapt colormap to the clustering. Default to \code{TRUE}.
#' @param ... additional parameters for S3 compatibility. Not used
#' @export
plot.ScoreSBM = function(x, type = c('data', 'predictedNetwork'), ordered = TRUE, color = TRUE, ...){
  x$plot(type, ordered, color)
}

#' Model Predictions
#'
#' Make predictions from a Score  SBM.
#'
#' @param object an R6 object inheriting from class ScoreSBM_fit
#' @param ... additional parameters for S3 compatibility. Not used
#' @return a matrix of expected values for each dyad of the underlying network G
#' @importFrom stats predict
#' @export
predict.ScoreSBM <- function(object, ...) {
  stopifnot(is_ScoreSBM(object))
  object$predict()
}


