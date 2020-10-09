#' R6 class for Bipartite SBM sampler
#'
#' @import R6
#' @include R6Class-SBM_sampler.R
#' @export
BipartiteSBM_sampler <-
  R6::R6Class(classname = "BipartiteSBM_sampler",
    inherit = SBM_sampler,
    public = list(
      #' @description constructor for SBM
      #' @param model character describing the type of model
      #' @param nbNodes number of nodes in the network
      #' @param blockProp parameters for block proportions (vector of list of vectors)
      #' @param connectParam list of parameters for connectivity with a matrix of means 'mean' and an optional scalar for the variance 'var'. The dimensions of mu must match \code{blockProp} lengths
      #' @param dimLabels optional labels of each dimension (in row, in column)
      #' @param covarParam optional vector of covariates effect
      #' @param covarList optional list of covariates data
      initialize = function(model, nbNodes, blockProp, connectParam, dimLabels=list(row="rowLabel", col="colLabel"), covarParam=numeric(0), covarList=list()) {
        ## SANITY CHECKS
        stopifnot(length(blockProp) ==  2,
                  length(blockProp[[1]]) ==  nrow(connectParam$mean),
                  length(blockProp[[2]]) ==  ncol(connectParam$mean))
        names(blockProp) <- c("row", "col")
        super$initialize(model, nbNodes, blockProp, dimLabels, connectParam, covarParam, covarList)
        self$rIncidence()
      },
      #' @description a method to generate a vector of block indicators
      rMemberships = function() {
        private$Z <- list(
          row = t(rmultinom(private$dim[1], size = 1, prob = private$pi[[1]])),
          col = t(rmultinom(private$dim[2], size = 1, prob = private$pi[[2]]))
          )
      },
      #' @description a method to sample an adjacency matrix for the current SBM
      #' @return nothing (sampled matrix is store in the current object, accessible via $netMatrix)
      rIncidence = function() {
        Y <- private$sampling_func(private$dim[1]*private$dim[2], list(mean = self$expectation, var = self$variance)) %>%
          matrix(private$dim[1], private$dim[2])
        private$Y <- Y
      },
      #' @description show method
      #' @param type character used to specify the type of SBM
      show = function(type = "Sampler for a Bipartite Stochastic Block Model") {
        super$show(type)
        cat("* R6 methods \n")
        cat("  $rMemberships(), $rIncidence() \n")
      }
    ),
    active = list(
      #' @field nbNodes vector of size 2: number of nodes (rows, columns)
      nbNodes     = function(value) {private$dim},
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
        if (self$nbCovariates > 0) mu <- private$invlink(private$link(mu) + self$covarEffect)
        mu
      }
    )
  )


