#' R6 Class definition of an Bipartite SBM fit
#'
#' This class is designed to give a representation and adjust an LBM fitted with blockmodels.
#'
#' @import R6 blockmodels
#' @include R6Class-SBM_fit.R
#' @export
BipartiteSBM_fit <-
  R6::R6Class(classname = "BipartiteSBM_fit",
    inherit = SBM_fit,
    public = list(
      #' @description constructor for a Bipartite SBM fit
      #' @param incidenceMatrix rectangular (weighted) matrix
      #' @param model character (\code{'bernoulli'}, \code{'poisson'}, \code{'gaussian'})
      #' @param covarList and optional list of covariates, each of whom must have the same dimension as \code{incidenceMatrix}
      initialize = function(incidenceMatrix, model, covarList=list()) {
        ## INITIALIZE THE SBM OBJECT ACCORDING TO THE DATA
        super$initialize(incidenceMatrix, model, covarList)

      },
      #' @description function to perform optimization
      #' @param verbosity integer, the level of verbosity. Default to 3
      #' @param plot logical, if TRUE ploting is done dynamically on the screen. Default to \code{TRUE}
      #' @param nbCores integer, the number of cores to use. Default is \code{parallel::detectCores()}.
      #' @param explorFactor double factor for exploring successive model
      #' @param nbBlocksRange 2-size vector: range of exploration
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
        args <- list(membership_type = "LBM", adj = private$Y)
        if (self$nbCovariates > 0) args$covariates <- private$X
        args <- c(args, blockmodelsOptions)

        ## model construction
        model_type <- ifelse(self$nbCovariates > 0, paste0(model,"_covariates"), private$model)
        BMobject <- do.call(paste0("BM_", model_type), args)

        ## performing estimation
        BMobject$estimate()

        ## Exporting blockmodels output to BipartiteSBM_fit fields
        ind_best      <- which.max(BMobject$ICL)
        private$J     <- BMobject$PL[ind_best]
        private$vICL  <- BMobject$ICL[ind_best]
        private$Y_hat <- BMobject$prediction(ind_best)
        private$tau   <- list(row = BMobject$memberships[[ind_best]]$Z1, col = BMobject$memberships[[ind_best]]$Z2)
        private$pi    <- lapply(private$tau, colMeans)
        parameters    <- BMobject$model_parameters[[ind_best]]
        private$beta  <- parameters$beta ## NULL if no covariates

        private$theta <- switch(model_type,
          "bernoulli"           = list(mean = parameters$pi),
          "bernoull_covariates" = list(mean = .logistic(parameters$m)),
          "poisson"             = list(mean = parameters$lambda),
          "poisson_covariates"  = list(mean = parameters$lambda),
          "gaussian"            = list(mean = parameters$mu, var = parameters$sigma2),
          "gaussian_covariates" = list(mean = parameters$mu, var = parameters$sigma2)
        )

        ## record fitted/expected value
        private$Y_hat <- self$predict()

        invisible(BMobject)
      },
      #' @description prediction under the currently estimated model
      #' @param covarList a list of covariates. By default, we use the covariates with which the model was estimated.
      predict = function(covarList = self$covarList) {
        stopifnot(is.list(covarList), self$nbCovariates == length(covarList))
        if (length(covarList) > 0) {
          stopifnot(all.equal(self$dimension[1], sapply(covarList, nrow)),
                    all.equal(self$dimension[2], sapply(covarList, ncol)))
        }
        mu <- private$tau[[1]] %*% private$theta$mean %*% t(private$tau[[2]])
        if (length(self$covList) > 0) mu <- private$invlink(private$link(mu) + self$covarEffect)
        mu
      },
      #' @description show method
      #' @param type character used to specify the type of SBM
      show = function(type = "Fit of a Bipartite Stochastic Block Model") super$show(type)
    ),
    active = list(
      #' @field nbNodes vector of size 2: number of nodes (rows, columns)
      nbNodes     = function(value) {private$dim},
      #' @field nbBlocks vector of size 2: number of blocks (rows, columns)
      nbBlocks    = function(value) {sapply(private$pi, length)},
      #' @field nbDyads number of dyads (potential edges in the network)
      nbDyads     = function(value) {private$dim[1] * private$dim[2]},
      #' @field memberships list of size 2: vector of memberships in row, in column.
      memberships = function(value) {lapply(private$tau, as_clustering)}
    )
  )
