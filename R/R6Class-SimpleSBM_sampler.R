#' R6 class for Simple SBM sampler
#'
#' @import R6
#' @include R6Class-SBM_sampler.R
#' @export
SimpleSBM_sampler <-
  R6::R6Class(classname = "SimpleSBM_sampler",
   inherit = SBM_sampler,
   ## fields for internal use (referring to the mathematical notation)
    private = list(
      directed_ = NULL # is the network directed or not (Symmetric netMatrix)
    ),
    public = list(
      #' @description constructor for SBM
      #' @param model character describing the type of model
      #' @param nbNodes number of nodes in the network
      #' @param directed logical, directed network or not.
      #' @param blockProp parameters for block proportions (vector of list of vectors)
      #' @param connectParam list of parameters for connectivity with a matrix of means 'mean' and an optional scalar for the variance 'var'. The size of mu must match \code{blockProp} length
      #' @param dimLabels optional labels of each dimension (in row, in column)
      #' @param covarParam optional vector of covariates effect
      #' @param covarList optional list of covariates data
      initialize = function(model, nbNodes, directed, blockProp, connectParam, dimLabels=list(row="rowLabel", col="colLabel"), covarParam=numeric(0), covarList=list()) {

        ## ADDITIONAL SANITY CHECKS
        stopifnot(all(blockProp > 0))
        stopifnot(all.equal(length(blockProp),      # dimensions match between vector of
                            ncol(connectParam$mean),  # block proportion and connectParam$mean
                            nrow(connectParam$mean)))
        if (!directed) stopifnot(isSymmetric(connectParam$mean)) # connectivity and direction must agree
        super$initialize(model, c(nbNodes, nbNodes), blockProp, dimLabels, connectParam, covarParam, covarList)
        private$directed_ <- directed
        self$rAdjacency()
      },
      #' @description a method to generate a vector of block indicators
      #' @return nothing (sampled memberships is stored in the current object)
      rMemberships = function() {
        private$Z <- t(rmultinom(private$dim[1], size = 1, prob = private$pi))
      },
      #' @description a method to sample an adjacency matrix for the current SBM
      #' @return nothing (sampled adjacency matrix is stored in the current object)
      rAdjacency = function() {
        Y <- suppressWarnings(private$sampling_func(self$nbNodes**2, list(mean = self$expectation, var = self$variance))) %>%
          matrix(private$dim[1], private$dim[2])
        diag(Y) <- NA
        if (!private$directed_) Y <- Y * lower.tri(Y) + t(Y * lower.tri(Y))
        private$Y <- Y
      },
      #' @description show method
      #' @param type character used to specify the type of SBM
      show = function(type = "Sampler for a Simple Stochastic Block Model") {
        super$show(type)
        cat("* R6 methods \n")
        cat("  $rMemberships(), $rAdjacency() \n")
      }
    ),
    active = list(
      #' @field nbNodes number of nodes
      nbNodes     = function(value) {private$dim[1]},
      #' @field nbBlocks number of blocks
      nbBlocks    = function(value) {length(private$pi)},
      #' @field nbDyads number of dyads (potential edges in the network)
      nbDyads     = function(value) {ifelse(private$directed_, self$nbNodes*(self$nbNodes - 1), self$nbNodes*(self$nbNodes - 1)/2)},
      #' @field nbConnectParam number of parameter used for the connectivity
      nbConnectParam = function(value) {ifelse(private$directed_, self$nbBlocks^2, self$nbBlocks*(self$nbBlocks + 1)/2)},
      #' @field memberships vector of clustering
      memberships = function(value) {if (!is.null(private$Z)) as_clustering(private$Z)},
      #' @field indMemberships matrix for clustering memberships
      indMemberships = function(value) {private$Z},
      #' @field expectation expected values of connection under the current model
      expectation = function() {
        mu <- private$Z %*% private$theta$mean %*% t(private$Z)
        if (self$nbCovariates > 0) mu <- private$invlink(private$link(mu) + self$covarEffect)
        diag(mu) <- NA
        mu
      },
      #' @field directed is the network directed or not
      directed = function(value) {private$directed_}
    )
  )

