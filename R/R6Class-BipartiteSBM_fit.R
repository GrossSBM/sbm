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
    private = list(
      import_from_BM  = function(index = which.max(private$BMobject$ICL)) {
        super$import_from_BM(index)
        private$tau <- list(
          row = private$BMobject$memberships[[index]]$Z1,
          col = private$BMobject$memberships[[index]]$Z2
        )
        private$pi  <- lapply(private$tau, colMeans)
      }
    ),
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
        args <- list(membership_type = "LBM", adj = private$Y)
        if (self$nbCovariates > 0) args$covariates <- private$X
        args <- c(args, blockmodelsOptions)

        ## model construction
        model_type <- ifelse(self$nbCovariates > 0, paste0(model,"_covariates"), private$model)
        if (model_type == 'bernoulli_covariates' & fast == TRUE) model_type <- 'bernoulli_covariates_fast'
        private$BMobject <- do.call(paste0("BM_", model_type), args)

        ## performing estimation
        private$BMobject$estimate()

        ## Exporting blockmodels output to BipartiteSBM_fit fields
        private$import_from_BM()

        invisible(private$BMobject)
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
      #' @description permute group labels by order of decreasing probability
      reorder = function() {
        oRow <- order(self$connectParam$mean %*% self$blockProp$col, decreasing = TRUE)
        oCol <- order(self$blockProp$row %*% self$connectParam$mean, decreasing = TRUE)
        private$pi[[1]] <- private$pi[[1]][oRow]
        private$pi[[2]] <- private$pi[[2]][oCol]
        private$theta$mean <- private$theta$mean[oRow, oCol]
        private$tau[[1]] <- private$tau[[1]][, oRow, drop = FALSE]
        private$tau[[2]] <- private$tau[[2]][, oCol, drop = FALSE]
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
      memberships = function(value) {lapply(private$tau, as_clustering)},
      #' @field storedModels data.frame of all models fitted (and stored) during the optimization
      storedModels = function(value) {
        rowBlocks <- c(0, unlist(sapply(private$BMobject$memberships, function(m) ncol(m$Z1))))
        colBlocks <- c(0, unlist(sapply(private$BMobject$memberships, function(m) ncol(m$Z2))))
        nbConnectParam <- c(NA, unlist(sapply(private$BMobject$model_parameters, function(param) param$n_parameters)))
        data.frame(
          nbParams  = nbConnectParam + rowBlocks + colBlocks - 2,
          rowBlocks = rowBlocks,
          colBlocks = colBlocks,
          nbBlocks  = rowBlocks + colBlocks,
          ICL       = private$BMobject$ICL,
          loglik    = private$BMobject$PL
        )
      }
    )
  )
