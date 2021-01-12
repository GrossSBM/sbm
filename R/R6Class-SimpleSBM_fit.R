#' R6 Class definition of a Simple SBM fit
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
      #--------------------------------------------
      #' @description constructor for a Simple SBM fit
      #' @param adjacencyMatrix square (weighted) matrix
      #' @param model character (\code{'bernoulli'}, \code{'poisson'}, \code{'gaussian'})
      #' @param directed logical, directed network or not. In not, \code{adjacencyMatrix} must be symmetric.
      #' @param dimLabels list of labels of each dimension (in row, in columns)
      #' @param covarList and optional list of covariates, each of whom must have the same dimension as \code{adjacencyMatrix}
      initialize = function(adjacencyMatrix, model, directed, dimLabels=list(row="node", col="node"), covarList=list()) {

        ## SANITY CHECKS
        stopifnot(all.equal(nrow(adjacencyMatrix), ncol(adjacencyMatrix)))  # matrix must be square
        stopifnot(isSymmetric(adjacencyMatrix) == !directed)                # symmetry and direction must agree

        ## INITIALIZE THE SBM OBJECT ACCORDING TO THE DATA
        super$initialize(adjacencyMatrix, model, dimLabels, covarList)
        private$directed_ <- directed

      },
      #--------------------------------------------
      #' @description function to perform optimization
      #' @param estimOptions a list of parameters controlling the inference algorithm and model selection. See details.
      #'
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

        if(private$model == 'ZIgaussian') stop("Inference not yet implemented for ZI gaussian network")

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
        args <- list(membership_type =  ifelse(!private$directed_, "SBM_sym", "SBM"), adj = .na2zero(private$Y))
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
      #--------------------------------------------
      #' @description prediction under the currently parameters
      #' @param covarList a list of covariates. By default, we use the covariates with which the model was estimated
      #' @param theta_p0 a threshold...
      #' @return a matrix of expected values for each dyad
      predict = function(covarList = self$covarList, theta_p0 = 0) {
        stopifnot(is.list(covarList), self$nbCovariates == length(covarList))
        mu <- private$tau %*% ( ((1-theta_p0)>0.5) * private$theta$mean ) %*% t(private$tau)
        if (self$nbCovariates > 0) {
          stopifnot(all(sapply(covarList, nrow) == self$nbNodes,
                        sapply(covarList, ncol) == self$nbNodes))
          mu <- private$invlink(private$link(mu) + self$covarEffect)
        }
        mu
      },
      #' @description permute group labels by order of decreasing probability
      reorder = function(){
        o <- order(private$theta$mean %*% private$pi, decreasing = TRUE)
        private$pi <- private$pi[o]
        private$theta$mean <- private$theta$mean[o,o]
        private$tau <- private$tau[, o, drop = FALSE]
      },
      #--------------------------------------------
      #' @description show method
      #' @param type character used to specify the type of SBM
      show = function(type = "Fit of a Simple Stochastic Block Model")
        {super$show(type)}
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
      #' @field blockProp vector of block proportions (aka prior probabilities of each block)
      blockProp   = function(value) {
        if (missing(value))
          return(private$pi)
        else {
          stopifnot(is.numeric(value))
          private$pi <- value
        }
      },
      #' @field probMemberships  matrix of estimated probabilities for block memberships for all nodes
      probMemberships = function(value) {
        if (missing(value))
          return(private$tau)
        else {
          stopifnot(nrow(value)==private$dim[1])
          private$tau <- value
        }
      },
      #' @field nbConnectParam number of parameter used for the connectivity
      nbConnectParam = function(value) {ifelse(private$directed_, self$nbBlocks^2, self$nbBlocks*(self$nbBlocks + 1)/2)},
      #' @field directed is the network directed or not
      directed = function(value) {private$directed_},
      #' @field penalty double, value of the penalty term in ICL
      penalty  = function(value) {(self$nbConnectParam + self$nbCovariates) * log(self$nbDyads) + (self$nbBlocks-1) * log(self$nbNodes)},
      #' @field entropy double, value of the entropy due to the clustering distribution
      entropy  = function(value) {-sum(.xlogx(private$tau))},
      #' @field storedModels data.frame of all models fitted (and stored) during the optimization
      storedModels = function(value) {
        nbBlocks <- unlist(sapply(private$BMobject$memberships, function(m) ncol(m$Z)))
        nbConnectParam <- unlist(sapply(private$BMobject$model_parameters, function(param) param$n_parameters))
        U <- data.frame(
          indexModel  = nbBlocks,
          nbParams = nbConnectParam + nbBlocks - 1,
          nbBlocks = nbBlocks,
          ICL      = private$BMobject$ICL,
          loglik   = private$BMobject$PL
          )
        U[!is.na(U$nbParams),]
      }
    )
  )
