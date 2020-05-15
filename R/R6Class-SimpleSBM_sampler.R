#' R6 class for Simple SBM sampler
#'
#' @import R6
#' @include R6Class-SBM.R
#' @export
SimpleSBM_sampler <- # this virtual class is the mother of all subtypes of SBM (Simple or Bipartite)
  R6::R6Class(classname = "SimpleSBM_sampler",
    inherit = SBM,
    ## fields for internal use (referring to the mathematical notation)
    private = list(
      directed_ = NULL, # is the network directed or not (Symmetric netMatrix)
      Z         = NULL, # the sampled indicator of blocks
      sampling_func = NULL #
    ),
    public = list(
      #' @description constructor for SBM
      #' @param model character describing the type of model
      #' @param nbNodes number of nodes in the network
      #' @param directed logical, directed network or not.
      #' @param blockProp parameters for block proportions (vector of list of vectors)
      #' @param connectParam list of parameters for connectivity with a matrix of means 'mu' and an optional matrix of variances 'sigma2', the sizes of which must match \code{blockProp} length
      #' @param covarParam optional vector of covariates effect
      #' @param covarList optional list of covariates data
      initialize = function(model, nbNodes, directed, blockProp, connectParam, covarParam=numeric(0), covarList=list()) {

        super$initialize(model = model, dimension = c(nbNodes, nbNodes), blockProp = blockProp, connectParam = connectParam, covarParam = covarParam, covarList = covarList)

        ## ADDITIONAL SANITY CHECKS
        stopifnot(all.equal(length(blockProp),      # dimensions match between vector of
                            ncol(connectParam$mu),  # block proportion and connectParam$mu
                            nrow(connectParam$mu)))
        stopifnot(isSymmetric(connectParam$mu) == !directed) # connectivity and direction must agree

        if (model == 'gaussian')
          stopifnot(all.equal(length(blockProp),      # dimensions match between vector of
                              ncol(connectParam$sigma2),  # block proportion and connectParam$mu
                              nrow(connectParam$sigma2)))

        private$sampling_func <- switch(model,
            "gaussian"  = function(n, param)  rnorm(n = n, mean   = param$mu, sd = sqrt(param$sigma2)) ,
            "poisson"   = function(n, param)  rpois(n = n, lambda = param$mu) ,
            "bernoulli" = function(n, param) rbinom(n = n, size = 1, prob   = param$mu),
          )
        private$directed_ <- directed

        self$rBlocks()
        self$rAdjMatrix()
      },
      #' @description a method to generate a vector of block indicators
      rBlocks = function() {
        private$Z <- t(rmultinom(private$dim[1], size = 1, prob = private$pi))
      },
      #' @description a method to sample an adjacency matrix for the current SBM
      #' @return nothing (sampled matrix is store in the current object, accessible via $netMatrix)
      rAdjMatrix = function() {
        Y <- private$sampling_func(self$nbNodes**2, list(mu = self$expectation, sigma2 = self$variance)) %>%
          matrix(private$dim[1], private$dim[2])
        if (!private$directed_) Y <- Y * lower.tri(Y) + t(Y * lower.tri(Y))
        private$Y <- Y
      }
    ),
    active = list(
      #' @field nbNodes number of nodes
      nbNodes     = function(value) {private$dim[1]},
      #' @field nbBlocks number of blocks
      nbBlocks    = function(value) {length(private$pi)},
      #' @field nbDyads number of dyads (potential edges in the network)
      nbDyads     = function(value) {ifelse(private$directed, self$nbNodes*(self$nbNodes - 1), self$nbNodes*(self$nbNodes - 1)/2)},
      #' @field memberships vector of clustering
      memberships = function(value) {if (!is.null(private$Z)) as_clustering(private$Z)},
      #' @field indMemberships matrix for clustering memberships
      indMemberships = function(value) {private$Z},
      #' @field expectation expected values of connection under the current model
      expectation = function() {
        mu <- private$Z %*% private$theta$mu %*% t(private$Z)
        if (self$nbCovariates > 0)
          mu <- private$invlink(private$link(mu) + self$covarEffect)
        mu
      },
      #' @field variance variances of each dyad under the current model
      variance = function() {
        if (private$model == 'gaussian')
          return(private$Z %*% private$theta$sigma2 %*% t(private$Z))
        else
          return(NULL)
      }
    )
  )

