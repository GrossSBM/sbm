#' R6 Class definition of an Bipartite SBM fit
#'
#' This class is designed to give a representation and adjust an LBM fitted with blockmodels.
#'
#' @include R6Class-SBM.R
#' @import R6 blockmodels
BipartiteSBM_fit <-
  R6::R6Class(classname = "BipartiteSBM_fit",
    inherit = SBM_fit,
    public = list(
      ## constructor
      initialize = function(incidenceMatrix, model, covarList=list()) {

        ### TODO - GET MORE CHECKS
        ## SANITY CHECKS
        stopifnot(all(sapply(covarList, nrow) == nrow(incidenceMatrix))) # all covariate matrices match the
        stopifnot(all(sapply(covarList, ncol) == ncol(incidenceMatrix))) # dimension of the adjancecy matrix

        ## INITIALIZE THE SBM OBJECT ACCORDING TO THE DATA
        super$initialize(data = incidenceMatrix, model = model, covarList = covarList)

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
        args <- list(membership_type = "LBM", adj = private$Y)
        if (self$nbCovariates > 0) args$covariates <- private$X
        args <- c(args, blockmodelsOptions)

        ## model construction
        BMobject <- do.call(private$optimizer_name, args)

        ## performing estimation
        BMobject$estimate()

        ## Exporting blockmodels output to BipartiteSBM_fit fields
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
      nbNodes         = function(value) {private$dim},
      nbBlocks        = function(value) {sapply(private$pi, length)},
      nbDyads         = function(value) {private$dim[1] * private$dim[2]},
      memberships     = function(value) {lapply(private$tau, as_clustering)}
    )
  )

