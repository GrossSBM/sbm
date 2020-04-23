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
      directed = NA
    ),
    public = list(
      ## constructor
      initialize = function(adjacencyMatrix, model, directed, covarList=list()) {

        ### TODO - GET MORE CHECKS
        ## SANITY CHECKS
        stopifnot(all.equal(nrow(adjacencyMatrix), ncol(adjacencyMatrix)))  # matrix must be square
        stopifnot(isSymmetric(adjacencyMatrix) == !directed)                # symmetry and direction must agree
        stopifnot(all(c(sapply(covarList, nrow),                            # all covariate matrices match the
                        sapply(covarList, ncol)) == nrow(adjacencyMatrix))) # dimension of the adjancecy matrix

        ## INITIALIZE THE SBM OBJECT ACCORDING TO THE DATA
        super$initialize(adjacencyMatrix, model, covarList)
        private$directed <- directed

      },
      ## optimizer = call to blockmodels
      optimize = function(verbosity     = 3,
                          plot          = "",
                          explorFactor  = 1.5,
                          nbBlocksRange = c(4,Inf),
                          nbCores       = parallel::detectCores()) {

        ## translate to blockmodels list of options
        blockmodelsOptions <- list(
          verbosity          = verbosity,
          plotting           = plot,
          explore_min        = nbBlocksRange[1],
          explore_max        = nbBlocksRange[2],
          ncores             = nbCores,
          exploration_factor = explorFactor
        )

        ## generating arguments for blockmodels call
        args <- list(membership_type =  ifelse(private$directed, "SBM_sym", "SBM"), adj = private$Y)
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
      ## active binding to access fields outside the class
      nbNodes     = function(value) {private$dim[1]},
      nbBlocks    = function(value) {length(private$pi)},
      nbDyads     = function(value) {ifelse(private$directed, self$nbNodes*(self$nbNodes - 1), self$nbNodes*(self$nbNodes - 1)/2)},
      memberships = function(value) {as_clustering(private$tau)}
    )
  )
