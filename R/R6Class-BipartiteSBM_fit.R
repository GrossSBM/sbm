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
      #' @param dimLabels labels of each dimension (in row, in columns)
      #' @param covarList and optional list of covariates, each of whom must have the same dimension as \code{incidenceMatrix}
      initialize = function(incidenceMatrix, model, dimLabels=list(row="rowLabel", col="colLabel"), covarList=list()) {
        ## INITIALIZE THE SBM OBJECT ACCORDING TO THE DATA
        super$initialize(incidenceMatrix, model, dimLabels, covarList)

      },
      #' @description function to perform optimization
      #' @param estimOptions a list of parameters controlling the inference algorithm and model selection. See details.
      #' @details The list of parameters \code{estimOptions} essentially tunes the optimization process and the variational EM algorithm, with the following parameters
      #'  \itemize{
      #'  \item{"nbCores"}{integer for number of cores used. Default is 2}
      #'  \item{"verbosity"}{integer for verbosity (0, 1). Default is 1}
      #'  \item{"plot"}{boolean, should the ICL by dynamically plotted or not. Default is TRUE}
      #'  \item{"exploreFactor"}{control the exploration of the number of groups}
      #'  \item{"nbBlocksRange"}{minimal and maximal number or blocks explored}
      #'  \item{"fast"}{logical: should approximation be used for Bernoulli model with covariates. Default to \code{TRUE}}
      #' }
      optimize = function(estimOptions = list()){

        currentOptions <- list(
          verbosity     = 3,
          plot          = TRUE,
          explorFactor  = 1.5,
          nbBlocksRange = c(4,Inf),
          nbCores       = 2,
          fast          = TRUE
        )
        currentOptions[names(estimOptions)] <- estimOptions

        ## Transform estimOptions to a suited for blockmodels list of options
        blockmodelsOptions <- list(
          verbosity          = currentOptions$verbosity,
          plotting           = if(currentOptions$plot) character(0) else "",
          explore_min        = currentOptions$nbBlocksRange[1],
          explore_max        = currentOptions$nbBlocksRange[2],
          ncores             = currentOptions$nbCores,
          exploration_factor = currentOptions$explorFactor
        )
        fast  <- currentOptions$fast

        ## generating arguments for blockmodels call
        if(private$model == 'ZIgaussian') stop("Inference not  yet  implemented for Bipartite ZI gaussian network")

        args <- list(membership_type = "LBM", adj = private$Y)
        if (self$nbCovariates > 0) args$covariates <- private$X
        args <- c(args, blockmodelsOptions)

        ## model construction

        model_type <- ifelse(self$nbCovariates > 0, paste0(private$model,"_covariates"), private$model)
        if (model_type == 'bernoulli_covariates' & fast == TRUE) model_type <- 'bernoulli_covariates_fast'
        private$BMobject <- do.call(paste0("BM_", model_type), args)

        ## performing estimation
        private$BMobject$estimate()

        ## Exporting blockmodels output to BipartiteSBM_fit fields
        private$import_from_BM()

        invisible(private$BMobject)
      },
      #' @description prediction under the current parameters
      #' @param covarList a list of covariates. By default, we use the covariates with which the model was estimated.
      predict = function(covarList = self$covarList) {
        mu <- predict_lbm(self$dimension,
                          self$nbCovariates,
                          private$link,
                          private$invlink,
                          private$tau,
                          private$theta$mean,
                          self$covarEffect,
                          covarList,
                          private$theta$p0)
        mu
      },
      #' @description permute group labels by order of decreasing probability
      reorder = function() {
        O <- order_lbm(private$theta$mean,private$pi)
        oRow <-O$row
        oCol <-O$col
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
      #' @field nbConnectParam number of parameter used for the connectivity
      nbConnectParam = function(value) {self$nbBlocks[1]*self$nbBlocks[2]},
      #' @field memberships list of size 2: vector of memberships in row, in column.
      memberships = function(value) {lapply(private$tau, as_clustering)},
      #' @field penalty double, value of the penalty term in ICL
      penalty  = function(value) {(self$nbConnectParam + self$nbCovariates) * log(self$nbDyads) + (self$nbBlocks[1]-1) * log(self$nbNodes[1]) + (self$nbBlocks[2]-1) * log(self$nbNodes[2])},
      #' @field entropy double, value of the entropy due to the clustering distribution
      entropy  = function(value) {-sum(.xlogx(private$tau[[1]]))-sum(.xlogx(private$tau[[2]]))},
      #' @field storedModels data.frame of all models fitted (and stored) during the optimization
      storedModels = function(value) {
        rowBlocks <- c(0, unlist(sapply(private$BMobject$memberships, function(m) ncol(m$Z1))))
        colBlocks <- c(0, unlist(sapply(private$BMobject$memberships, function(m) ncol(m$Z2))))
        nbConnectParam <- c(NA, unlist(sapply(private$BMobject$model_parameters, function(param) param$n_parameters)))
        U <- data.frame(
          indexModel = rowBlocks + colBlocks,
          nbParams  = nbConnectParam + rowBlocks + colBlocks - 2,
          rowBlocks = rowBlocks,
          colBlocks = colBlocks,
          nbBlocks  = rowBlocks + colBlocks,
          ICL       = private$BMobject$ICL,
          loglik    = private$BMobject$PL
        )
        U[!is.na(U$nbParams),];

      }
    )
  )
