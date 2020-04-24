#' R6 Class definition of an Simple SBM fit
#'
#' This class is designed to give a representation and adjust an SBM fitted with blockmodels.
#'
#' @import R6 blockmodels
#' @include R6Class-SBM.R
SimpleSBM_fit <-
  R6::R6Class(classname = "SimpleSBM_fit",
    inherit = SBM_fit,
    private = list(
      directed_ = NULL
    ),
    public = list(
      ## constructor
      initialize = function(adjacencyMatrix, model, directed, covarList=list()) {

        ## SANITY CHECKS
        stopifnot(all.equal(nrow(adjacencyMatrix), ncol(adjacencyMatrix)))  # matrix must be square
        stopifnot(isSymmetric(adjacencyMatrix) == !directed)                # symmetry and direction must agree
        stopifnot(all(c(sapply(covarList, nrow),                            # all covariate matrices match the
                        sapply(covarList, ncol)) == nrow(adjacencyMatrix))) # dimension of the adjancecy matrix

        ## INITIALIZE THE SBM OBJECT ACCORDING TO THE DATA
        super$initialize(adjacencyMatrix, model, covarList)
        private$directed_ <- directed

      },
      #' @description function to perform optimization
      #' @param verbosity integer, the level of verbosity. Default to 3
      #' @param plot logical, if TRUE ploting is done dynamically on the screen. Default to \code{TRUE}
      #' @param nbCores integer, the number of cores to use. Default is \code{parallel::detectCores()}.
      optimize = function(verbosity     = 3,
                          plot          = TRUE,
                          explorFactor  = 1.5,
                          nbBlocksRange = c(4,Inf),
                          nbCores       = parallel::detectCores()) {

        ## translate to blockmodels list of options
        blockmodelsOptions <- list(
          verbosity          = verbosity,
          plotting           = ifelse(plot, character(0), ""),
          explore_min        = nbBlocksRange[1],
          explore_max        = nbBlocksRange[2],
          ncores             = nbCores,
          exploration_factor = explorFactor
        )

        ## generating arguments for blockmodels call
        args <- list(membership_type =  ifelse(private$directed_, "SBM_sym", "SBM"), adj = private$Y)
        if (self$nbCovariates > 0) args$covariates <- private$X
        args <- c(args, blockmodelsOptions)

        ## model construction
        BMobject <- do.call(private$optimizer_name, args)

        ## performing estimation
        BMobject$estimate()

        ## Exporting blockmodels output to simpleSBM_fit fields
        ind_best      <- which.max(BMobject$ICL)
        private$J     <- BMobject$PL[ind_best]
        private$vICL  <- BMobject$ICL[ind_best]
        private$Y_hat <- BMobject$prediction(ind_best)

        res <- extractParamBM(BMobject, ind_best)

        private$pi    <- res$blockProp
        private$theta <- res$connectParam
        private$beta  <- res$covarParam
        private$tau   <- res$probMemberships

        invisible(BMobject)
      }
    ),
    active = list(
      #' @field number of nodes
      nbNodes     = function(value) {private$dim[1]},
      #' @field number of blocks
      nbBlocks    = function(value) {length(private$pi)},
      #' @field number of dyads (potential edges in the network)
      nbDyads     = function(value) {ifelse(private$directed, self$nbNodes*(self$nbNodes - 1), self$nbNodes*(self$nbNodes - 1)/2)},
      #' @field vector of clustering
      memberships = function(value) {as_clustering(private$tau)}
    )
  )

#' R6 class for Simple SBM sampler
#'
#' @import R6
SimpleSBM_sampler <- # this virtual class is the mother of all subtypes of SBM (Simple or Bipartite)
  R6::R6Class(classname = "SimpleSBM_sampler",
    inherit = SBM,
    ## fields for internal use (referring to the mathematical notation)
    private = list(
      Z            = NULL  # the sampled indicator of blocks
    ),
    public = list(
      #' @description a method to generate a vector of clusters indicators
      rBlocks = function() {
        private$Z <- t(rmultinom(private$dim[1], size = 1, prob = private$pi))
      }
      # ,
      # ## a method to sample an adjacency matrix for the current SBM
      # rAdjMatrix = function() {
      #
      #
      #   Y <- matrix(rbinom(private$N^2, 1, self$connectProb), private$N)
      #
      #   if (!private$directed) Y <- Y * lower.tri(Y) + t(Y * lower.tri(Y))
      #
      #   private$Y <- Y
      # }
    ),
    active = list(
      indMemberships = function(value) {private$Z}
      # ,
      # connectProb = function(value) {
      #   PI <- private$Z %*% private$theta %*% t(private$Z)
      #   if (self$nbCovariates > 0) {
      #     PI <- logistic(PI + roundProduct(simplify2array(private$X), private$beta))
      #   }
      #   PI
      # }
    )
  )

