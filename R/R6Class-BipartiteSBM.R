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

        ### TODO - GET MORE CHECKS
        ## SANITY CHECKS
        stopifnot(all(sapply(covarList, nrow) == nrow(incidenceMatrix))) # all covariate matrices match the
        stopifnot(all(sapply(covarList, ncol) == ncol(incidenceMatrix))) # dimension of the adjancecy matrix

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
        BMobject <- do.call(paste0("BM_", private$model), args)

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

        private$theta <- switch(private$model,
          "bernoulli"           = list(mu = parameters$pi),
          "bernoull_covariates" = list(mu = .logistic(parameters$m)),
          "poisson"             = list(mu = parameters$lambda),
          "poisson_covariates"  = list(mu = parameters$lambda),
          "gaussian"            = list(mu = parameters$mu, sigma = parameters$sigma2),
          "gaussian_covariates" = list(mu = parameters$mu, sigma = parameters$sigma2)
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
        mu <- private$tau[[1]] %*% private$theta$mu %*% t(private$tau[[2]])
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

#' R6 class for Simple SBM sampler
#'
#' @import R6
BipartiteSBM_sampler <-
  R6::R6Class(classname = "BipartiteSBM_sampler",
    inherit = SBM,
    ## fields for internal use (referring to the mathematical notation)
    private = list(
      Z            = NULL  # the sampled indicator of blocks
    ),
    public = list(
      #' @description a method to generate a vector of clusters indicators
      rBlocks = function() {
        private$Z <- list(
          row = t(rmultinom(private$dim[1], size = 1, prob = private$pi[[1]])),
          col = t(rmultinom(private$dim[2], size = 1, prob = private$pi[[2]]))
          )

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
    )#,
    # active = list(
    #   indMemberships = function(value) {private$Z}
    #   # ,
    # )
  )

