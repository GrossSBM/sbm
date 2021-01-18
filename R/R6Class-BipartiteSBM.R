#' R6 class for Bipartite SBM
#'
#' @import R6
#' @include R6Class-SBM.R
#' @export
BipartiteSBM <-
  R6::R6Class(classname = "BipartiteSBM",
    inherit = SBM,
    public = list(
      #' @description constructor for SBM
      #' @param model character describing the type of model
      #' @param nbNodes number of nodes in the network
      #' @param blockProp parameters for block proportions (vector of list of vectors)
      #' @param connectParam list of parameters for connectivity with a matrix of means 'mean' and an optional scalar for the variance 'var'. The dimensions of mu must match \code{blockProp} lengths
      #' @param dimLabels optional labels of each dimension (in row, in column)
      #' @param covarParam optional vector of covariates effect
      #' @param covarList optional list of covariates data
      initialize = function(model, nbNodes, blockProp, connectParam, dimLabels=c(row="row", col="col"), covarParam=numeric(0), covarList=list()) {

        ## SANITY CHECKS
        stopifnot(length(dimLabels) == 2)
        stopifnot(length(blockProp) ==  2, is.list(blockProp),
                  length(blockProp[[1]]) ==  nrow(connectParam$mean), # dimensions match between vector of
                  length(blockProp[[2]]) ==  ncol(connectParam$mean)) # block proportion and connectParam$mean
        stopifnot(all(blockProp[[1]] > 0), all(blockProp[[1]] < 1))   # positive proportions
        stopifnot(all(blockProp[[2]] > 0), all(blockProp[[2]] < 1))
        names(blockProp) <- names(dimLabels)
        super$initialize(model, NA, nbNodes, dimLabels, blockProp, connectParam, covarParam, covarList)
      },
      #' @description a method to generate a vector of block indicators
      #' @param store should the sampled blocks be stored (and overwrite the existing data)? Default to FALSE
      #' @return the sampled blocks
      rMemberships = function(store = FALSE) {
        Z <- list(
          row = t(rmultinom(private$dim[1], size = 1, prob = private$pi[[1]])),
          col = t(rmultinom(private$dim[2], size = 1, prob = private$pi[[2]]))
          )
        if (store) private$Z <- Z
        Z
      },
      #' @description a method to sample a network data for the current SBM
      #' @param store should the sampled network be stored (and overwrite the existing data)? Default to FALSE
      #' @return the sampled network
      rIncidence = function(store = FALSE) {
        Y <- private$sampling_func[[1]](private$dim[1]*private$dim[2], list(mean = self$expectation, var = private$theta$var)) %>%
          matrix(private$dim[1], private$dim[2])
        if (store) private$Y <- Y
        Y
      },
      #' @description show method
      #' @param type character used to specify the type of SBM
      show = function(type = "Bipartite Stochastic Block Model") {
        super$show(type)
        cat("* R6 methods \n")
        cat("  $rMemberships(), $rIncidence() \n")
        cat("  $indMemberships, $memberships, $expectation\n")
        cat("* S3 methods \n")
        cat("  plot, print, coef \n")      }
    ),
    active = list(
      #' @field nbBlocks vector of size 2: number of blocks (rows, columns)
      nbBlocks    = function(value) {sapply(private$pi, length)},
      #' @field nbDyads number of dyads (potential edges in the network)
      nbDyads     = function(value) {private$dim[1] * private$dim[2]},
      #' @field nbConnectParam number of parameter used for the connectivity
      nbConnectParam = function(value) {self$nbBlocks^2},
      #' @field memberships list of size 2: vector of memberships in row, in column.
      memberships = function(value) {lapply(private$Z, as_clustering)},
      #' @field indMemberships list of 2 matrix for clustering memberships
      indMemberships = function(value) {private$Z},
      #' @field expectation expected values of connection under the current model
      expectation = function() {
        mu <- private$Z[[1]] %*% private$theta$mean %*% t(private$Z[[2]])
        if (self$nbCovariates > 0) mu <- private$invlink[[1L]](private$link[[1L]](mu) + self$covarEffect)
        mu
      },
      #' @field variance variance of each dyad under the current model
      variance = function() {if (private$model == 'gaussian') return(private$theta$var) else return(NULL) }
    )
  )


