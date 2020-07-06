#' R6 Class definition of an Simple SBM fit
#'
#' This class is designed to give a representation and adjust an SBM fitted with blockmodels.
#'
#' @import R6 blockmodels
#' @include R6Class-SBM_fit.R
#' @export
SimpleSBM_fit <-
  R6::R6Class(classname = "SimpleSBM_fit",
    inherit = SBM_fit,
    private = list(
      directed_ = NULL, # is the network directed or not
      import_from_BM = function(index = which.max(private$BMobject$ICL)) { # a function updating the Class
        super$import_from_BM(index)
        private$tau <- private$BMobject$memberships[[index]]$Z
        private$pi  <- colMeans(private$tau)
      }
    ),
    public = list(
      #' @description constructor for a Simple SBM fit
      #' @param adjacencyMatrix square (weighted) matrix
      #' @param model character (\code{'bernoulli'}, \code{'poisson'}, \code{'gaussian'})
      #' @param directed logical, directed network or not. In not, \code{adjacencyMatrix} must be symmetric.
      #' @param covarList and optional list of covariates, each of whom must have the same dimension as \code{adjacencyMatrix}
      initialize = function(adjacencyMatrix, model, directed, covarList=list()) {

        ## SANITY CHECKS
        stopifnot(all.equal(nrow(adjacencyMatrix), ncol(adjacencyMatrix)))  # matrix must be square
        stopifnot(isSymmetric(adjacencyMatrix) == !directed)                # symmetry and direction must agree

        ## INITIALIZE THE SBM OBJECT ACCORDING TO THE DATA
        super$initialize(adjacencyMatrix, model, covarList)
        private$directed_ <- directed

      },
      #' @description function to perform optimization
      #' @param verbosity integer, the level of verbosity. Default to 3
      #' @param plot logical, if TRUE ploting is done dynamically on the screen. Default to \code{TRUE}
      #' @param nbCores integer, the number of cores to use. Default is 2.
      #' @param explorFactor double factor for exploring successive model
      #' @param nbBlocksRange 2-size vector: range of exploration
      #' @param fast logical: should approximation be used for Bernoulli model with covariates. Default to \code{TRUE}
      optimize = function(verbosity     = 3,
                          plot          = FALSE,
                          explorFactor  = 1.5,
                          nbBlocksRange = c(4,Inf),
                          nbCores       = 2,
                          fast          = TRUE) {

        ## translate to blockmodels list of options
        blockmodelsOptions <- list(
          verbosity          = verbosity,
          plotting           = if(plot) character(0) else "",
          explore_min        = nbBlocksRange[1],
          explore_max        = nbBlocksRange[2],
          ncores             = nbCores,
          exploration_factor = explorFactor
        )

        ## generating arguments for blockmodels call
        args <- list(membership_type =  ifelse(private$directed_, "SBM_sym", "SBM"), adj = .na2zero(private$Y))
        if (self$nbCovariates > 0) args$covariates <- private$X
        args <- c(args, blockmodelsOptions)

        ## model construction
        model_type <- ifelse(self$nbCovariates > 0, paste0(private$model,"_covariates"), private$model)
        if (model_type == 'bernoulli_covariates' & fast == TRUE) model_type <- 'bernoulli_covariates_fast'
        private$BMobject <- do.call(paste0("BM_", model_type), args)

        ## performing estimation
        private$BMobject$estimate()

        ## Exporting blockmodels output to simpleSBM_fit fields
        private$import_from_BM()

        invisible(private$BMobject)
      },
      #' @description prediction under the currently estimated model
      #' @param covarList a list of covariates. By default, we use the covariates with which the model was estimated
      #' @return a matrix of expected values for each dyad
      predict = function(covarList = self$covarList) {
        stopifnot(is.list(covarList), self$nbCovariates == length(covarList))
        mu <- private$tau %*% private$theta$mean %*% t(private$tau)
        if (self$nbCovariates > 0) {
          all(sapply(covarList, nrow) == self$nbNodes, sapply(covarList, ncol) == self$nbNodes)
          mu <- private$invlink(private$link(mu) + self$covarEffect)
        }
        mu
      },
      #' @description permute group labels by order of decreasing probability
      reorder = function() {
        o <- order(self$connectParam$mean %*% self$blockProp, decreasing = TRUE)
        private$pi <- private$pi[o]
        private$theta$mean <- private$theta$mean[o,o]
        private$tau <- private$tau[, o, drop = FALSE]
      },
      #' @description show method
      #' @param type character used to specify the type of SBM
      show = function(type = "Fit of a Simple Stochastic Block Model") super$show(type)
    ),
    active = list(
      #' @field nbNodes number of nodes
      nbNodes     = function(value) {private$dim[1]},
      #' @field nbBlocks number of blocks
      nbBlocks    = function(value) {length(private$pi)},
      #' @field nbDyads number of dyads (potential edges in the network)
      nbDyads     = function(value) {ifelse(private$directed_, self$nbNodes*(self$nbNodes - 1), self$nbNodes*(self$nbNodes - 1)/2)},
      #' @field memberships vector of clustering
      memberships = function(value) {as_clustering(private$tau)},
      #' @field directed is the network directed or not
      directed = function(value) {private$directed_},
      #' @field storedModels data.frame of all models fitted (and stored) during the optimization
      storedModels = function(value) {
        nbBlocks <- unlist(sapply(private$BMobject$memberships, function(m) ncol(m$Z)))
        nbConnectParam <- unlist(sapply(private$BMobject$model_parameters, function(param) param$n_parameters))
        data.frame(
          nbParams = nbConnectParam + nbBlocks - 1,
          nbBlocks = nbBlocks,
          ICL      = private$BMobject$ICL,
          loglik   = private$BMobject$PL
          )
      }
    )
  )
