#' R6 Class definition of an Bipartite SBM fit
#'
#' This class is designed to give a representation and adjust an LBM fitted with blockmodels.
#'
#' @import R6 blockmodels
#' @export
BipartiteSBM_fit <-
  R6::R6Class(
    classname = "BipartiteSBM_fit",
    inherit = BipartiteSBM,
    private = list(
      J              = NULL, # approximation of the log-likelihood
      vICL           = NULL, # approximation of the ICL
      tau            = NULL, # parameters for posterior probability of class belonging
      BMobject       = NULL, # blockmodels output (used to stored the optimization results when blockmodels is used)
      import_from_BM  = function(index = which.max(private$BMobject$ICL)) {
        private$J     <- private$BMobject$PL[index]
        private$vICL  <- private$BMobject$ICL[index]
        parameters    <- private$BMobject$model_parameters[[index]]
        private$beta  <- parameters$beta ## NULL if no covariates
        private$theta <- switch(private$BMobject$model_name,
          "bernoulli"                 = list(mean = parameters$pi),
          "bernoulli_covariates"      = list(mean = .logistic(parameters$m)),
          "bernoulli_covariates_fast" = list(mean = .logistic(parameters$m)),
          "poisson"                   = list(mean = parameters$lambda),
          "poisson_covariates"        = list(mean = parameters$lambda),
          "gaussian"                  = list(mean = parameters$mu, var = parameters$sigma2),
          "gaussian_covariates"       = list(mean = parameters$mu, var = parameters$sigma2),
          "ZIgaussian"                = list(mean = parameters$mu, var = parameters$sigma2, p0 = parameters$p0),
        )
        private$tau <- list(
          row = private$BMobject$memberships[[index]]$Z1,
          col = private$BMobject$memberships[[index]]$Z2
        )
        private$Z  <- private$tau
        private$pi <- lapply(private$tau, colMeans)
      }
    ),
    public = list(
      #' @description constructor for a Bipartite SBM fit
      #' @param incidenceMatrix rectangular (weighted) matrix
      #' @param model character (\code{'bernoulli'}, \code{'poisson'}, \code{'gaussian'})
      #' @param dimLabels labels of each dimension (in row, in columns)
      #' @param covarList and optional list of covariates, each of whom must have the same dimension as \code{incidenceMatrix}
      initialize = function(incidenceMatrix, model, dimLabels=c(row="rowName", col="colName"), covarList=list()) {

        ## SANITY CHECKS on data
        stopifnot(is.matrix(incidenceMatrix))                            # must be a matrix
        stopifnot(all(sapply(covarList, nrow) == nrow(incidenceMatrix))) # consistency of the covariates
        stopifnot(all(sapply(covarList, ncol) == ncol(incidenceMatrix))) # with the network data

        ## INITIALIZE THE SBM OBJECT ACCORDING TO THE DATA
        connectParam <- switch(model,
          "bernoulli"  = list(mean = matrix(0, 0, 0)),
          "poisson"    = list(mean = matrix(0, 0, 0)),
          "gaussian"   = list(mean = matrix(0, 0, 0), var = 1),
          "ZIgaussian" = list(mean = matrix(0, 0, 0), var = 1, p0 = 0),
        )

        ## INITIALIZE THE SBM OBJECT ACCORDING TO THE DATA
        super$initialize(model        = model,
                         nbNodes      = dim(incidenceMatrix),
                         blockProp    = rep(list(vector("numeric", 0)), 2),
                         connectParam = connectParam,
                         dimLabels    = dimLabels,
                         covarList    = covarList)
        private$Y <- incidenceMatrix
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

        if(private$model == 'ZIgaussian') stop("Inference not  yet  implemented for Bipartite ZI gaussian network")

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
        fast <- currentOptions$fast

        ## generating arguments for blockmodels call

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
      #' @param theta_p0 double for thresholding...
      predict = function(covarList = self$covarList, theta_p0 = 0) {
        stopifnot(!is.null(private$tau[[1]]),
                  !is.null(private$tau[[2]]),
                  !is.null(private$theta$mean))
        stopifnot(is.list(covarList),  self$nbCovariates == length(covarList))

        mu <- private$tau[[1]] %*% ( ((1-theta_p0)>0.5 ) * private$theta$mean )  %*% t(private$tau[[2]])
        if (length(covarList) > 0) {
          stopifnot(all(sapply(covarList, nrow) == self$dimension[1]),
                    all(sapply(covarList, ncol) == self$dimension[2]))
          mu <- private$invlink[[1L]](private$link[[1L]](mu) + self$covarEffect)
        }
        mu
      },
      #' @description method to select a specific model among the ones fitted during the optimization.
      #'  Fields of the current SBM_fit will be updated accordingly.
      #' @param index integer, the index of the model to be selected (row number in storedModels)
      setModel = function(index) {
        stopifnot(!is.null(private$BMobject))
        stopifnot(index %in% seq.int(nrow(self$storedModels)))
        private$import_from_BM(index)
        self$reorder()
      },
      #' @description permute group labels by order of decreasing probability
      reorder = function() {
        oRow <- order(private$theta$mean %*% private$pi[[2]], decreasing = TRUE)
        oCol <- order(private$pi[[1]] %*% private$theta$mean, decreasing = TRUE)
        private$pi[[1]] <- private$pi[[1]][oRow]
        private$pi[[2]] <- private$pi[[2]][oCol]
        private$theta$mean <- private$theta$mean[oRow, oCol]
        private$tau[[1]] <- private$tau[[1]][, oRow, drop = FALSE]
        private$tau[[2]] <- private$tau[[2]][, oCol, drop = FALSE]
      }
    ),
    active = list(
      #' @field memberships list of size 2: vector of memberships in row, in column.
      memberships = function(value) {lapply(private$tau, as_clustering)},
      #' @field blockProp list of block proportions (aka prior probabilities of each block)
      blockProp   = function(value) {
        if (missing(value))
          return(private$pi)
        else {
          stopifnot(is.list(value), length(value) == 2)
          private$pi <- value
        }
      },
      #' @field probMemberships list of 2 matrices of estimated probabilities for block memberships for all nodes
      probMemberships = function(value) {
        if (missing(value))
          return(private$tau)
        else {
          stopifnot(is.list(value))
          stopifnot(nrow(value[[1]]) == private$dim[1])
          stopifnot(nrow(value[[2]]) == private$dim[2])
          private$tau <- value
        }
      },
      #' @field connectParam parameters associated to the connectivity of the SBM, e.g. matrix of inter/inter block probabilities when model is Bernoulli
      connectParam   = function(value) {
        if (missing(value))
          return(private$theta)
        else {
          stopifnot(is.list(value))
          private$theta <- value
        }
      },
      #' @field loglik double: approximation of the log-likelihood (variational lower bound) reached
      loglik = function(value) {private$J},
      #' @field ICL double: value of the integrated classification log-likelihood
      ICL    = function(value) {private$vICL},
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
        U[!is.na(U$nbParams), , drop = FALSE]
      }
    )
  )
