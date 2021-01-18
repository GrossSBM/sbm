#' R6 class for Simple SBM
#'
#' @import R6
#' @include R6Class-SBM.R
#' @export
SimpleSBM <-
  R6::R6Class(classname = "SimpleSBM",
    inherit = SBM,
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
      initialize = function(model, nbNodes, directed, blockProp, connectParam, dimLabels=c(node="nodeName"), covarParam=numeric(0), covarList=list()) {

        ## ADDITIONAL SANITY CHECKS
        stopifnot(length(dimLabels) == 1)
        stopifnot(is.atomic(blockProp), all(blockProp > 0), all(blockProp < 1)) # positive proportions
        stopifnot(all.equal(length(blockProp),            # dimensions match between vector of
                            ncol(connectParam$mean),      # block proportion and connectParam$mean
                            nrow(connectParam$mean)))

        if (!directed) stopifnot(isSymmetric(connectParam$mean)) # connectivity and direction must agree
        super$initialize(model, directed, nbNodes, dimLabels, blockProp, connectParam, covarParam, covarList)
      },
      #' @description a method to generate a vector of block indicators
      #' @param store should the sampled blocks be stored (and overwrite the existing data)? Default to FALSE
      #' @return the sampled blocks
      rMemberships = function(store = FALSE) {
        Z <- t(rmultinom(private$dim, size = 1, prob = private$pi))
        if (store) private$Z <- Z
        Z
      },
      #' @description a method to sample a network data for the current SBM
      #' @param store should the sampled network be stored (and overwrite the existing data)? Default to FALSE
      #' @return the sampled network
      rAdjacency = function(store = FALSE) {
        Y <- suppressWarnings(private$sampling_func[[1]](self$nbNodes**2, list(mean = self$expectation, var = private$theta$var))) %>%
          matrix(private$dim, private$dim)
        diag(Y) <- NA
        if (!private$directed_) Y <- Y * lower.tri(Y) + t(Y * lower.tri(Y))
        if (store) private$Y <- Y
        Y
      },
      #' @description show method
      #' @param type character used to specify the type of SBM
      show = function(type = "Sampler for a Simple Stochastic Block Model") {
        super$show(type)
        cat("* R6 methods \n")
        cat("  $rMemberships(), $rAdjacency() \n")
        cat("  $indMemberships, $memberships, $expectation\n")
        cat("* S3 methods \n")
        cat("  plot, print, coef \n")
      }
    ),
    active = list(
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
        if (self$nbCovariates > 0) mu <- private$invlink[[1L]](private$link[[1L]](mu) + self$covarEffect)
        diag(mu) <- NA
        mu
      },
      #' @field variance variance of each dyad under the current model
      variance = function() {if (private$model == 'gaussian') return(private$theta$var) else return(NULL) }
    )
  )

