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
      #' @param dimLabels optional label for the node (default is "nodeName")
      #' @param covarParam optional vector of covariates effect
      #' @param covarList optional list of covariates data
      initialize = function(model, nbNodes, directed, blockProp, connectParam, dimLabels=c(node="nodeName"), covarParam=numeric(length(covarList)), covarList=list()) {

        ## SANITY CHECKS (on parameters)
        stopifnot(length(dimLabels) == 1)
        stopifnot(is.atomic(blockProp), all(blockProp > 0), all(blockProp < 1)) # positive proportions
        stopifnot(all.equal(length(blockProp), ncol(connectParam$mean)),        # dimensions match between vector of
                  all.equal(length(blockProp), nrow(connectParam$mean)))        # block proportion and connectParam$mean

        ## Check that connectivity parameters and model are consistent
        switch(model,
          "bernoulli"  = stopifnot(all(connectParam$mean >= 0), all(connectParam$mean <= 1)),
          "poisson"    = stopifnot(all(connectParam$mean >= 0)),
          "gaussian"   = stopifnot(length(connectParam$var) == 1, connectParam$var > 0),
          "ZIgaussian" = stopifnot(all(connectParam$p0 >= 0), all(connectParam$p0 <= 1))
        )

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
      #--------------------------------------------
      #' @description prediction under the currently parameters
      #' @param covarList a list of covariates. By default, we use the covariates with which the model was estimated
      #' @param theta_p0 a threshold...
      #' @return a matrix of expected values for each dyad
      predict = function(covarList = self$covarList, theta_p0 = 0) {
        stopifnot(is.list(covarList), self$nbCovariates == length(covarList))
        mu <- private$Z %*% ( ((1-theta_p0)>0.5) * private$theta$mean ) %*% t(private$Z)
        if (self$nbCovariates > 0) {
          stopifnot(all(sapply(covarList, nrow) == self$nbNodes,
                        sapply(covarList, ncol) == self$nbNodes))
          mu <- private$invlink[[1L]](private$link[[1L]](mu) + self$covarEffect)
        }
        mu
      },
      #' @description show method
      #' @param type character used to specify the type of SBM
      show = function(type = "Simple Stochastic Block Model") {
        super$show(type)
        cat("* R6 methods \n")
        cat("  $rMemberships(), $rAdjacency() \n")
        cat("  $indMemberships, $memberships, $expectation\n")
        cat("* S3 methods \n")
        cat("  plot, print, coef \n")
      },
      #' @description basic matrix plot method for SimpleSBM object or mesoscopic plot
      #' @param type character for the type of plot: either 'data' (true connection), 'expected' (fitted connection) or 'meso' (mesoscopic view). Default to 'data'.
      #' @param ordered logical: should the rows and columns be reordered according to the clustering? Default to \code{TRUE}.
      #' @param plotOptions list with the parameters for the plot. See help of the corresponding S3 method for details.
      #' @return a ggplot2 object for the \code{'data'} and \code{'expected'}, a list with the igraph object \code{g}, the \code{layout} and the \code{plotOptions} for the \code{'meso'}
      #' @import ggplot2
      plot = function(type = c('data','expected','meso'), ordered = TRUE, plotOptions = list()) {

        if(is.null(self$memberships)) {ordered <- FALSE; type <- 'data'}
        if (ordered & !is.null(self$memberships))
          clustering <- list(row = self$memberships)
        else
          clustering <- NULL

        switch(match.arg(type),
          "meso" =
            plotMeso(
              thetaMean  = private$theta$mean,
              pi         = private$pi,
              model      = private$model,
              directed   = private$directed_,
              bipartite  = FALSE,
              nbNodes    = self$dimension,
              nodeLabels = as.list(private$dimlab),
              plotOptions),
          "data" =
            plotMatrix(self$networkData,
                       private$dimlab,
                       clustering, plotOptions),
          "expected" =
            plotMatrix(self$expectation,
                       private$dimlab,
                       clustering, plotOptions)
        )
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
      indMemberships = function(value) {as_indicator(as_clustering(private$Z))}
    )
  )

