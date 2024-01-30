#' R6 Class definition of a Simple SBM fit
#'
#' This class is designed to give a representation and adjust an SBM fitted with blockmodels.
#'
#' @import R6 blockmodels
#' @export
SimpleSBM_fit <-
  R6::R6Class(
    classname = "SimpleSBM_fit",
    inherit = SimpleSBM,
    private = list(
      J              = NULL, # approximation of the log-likelihood
      vICL           = NULL, # approximation of the ICL
      BMobject       = NULL, # blockmodels output (used to stored the optimization results when blockmodels is used)
      import_from_BM = function(index = which.max(private$BMobject$ICL)) { # a function updating the Class
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
        private$Z  <- private$BMobject$memberships[[index]]$Z
        private$pi <- colMeans(private$Z)
      }
    ),
    public = list(
      #' @description constructor for a Simple SBM fit
      #' @param adjacencyMatrix square (weighted) matrix
      #' @param model character (\code{'bernoulli'}, \code{'poisson'}, \code{'gaussian'})
      #' @param directed logical, directed network or not. In not, \code{adjacencyMatrix} must be symmetric.
      #' @param dimLabels list of labels of each dimension (in row, in columns)
      #' @param covarList and optional list of covariates, each of whom must have the same dimension as \code{adjacencyMatrix}
      initialize = function(adjacencyMatrix, model, directed, dimLabels=c(node="nodeName"), covarList=list()) {

        ## SANITY CHECKS (on data)
        stopifnot(is.matrix(adjacencyMatrix))                   # must be a matrix
        stopifnot(all.equal(nrow(adjacencyMatrix),
                            ncol(adjacencyMatrix)))             # matrix must be square
        stopifnot(isSymmetric(adjacencyMatrix) == !directed)    # symmetry and direction must agree
        stopifnot(all(sapply(covarList, nrow) == nrow(adjacencyMatrix))) # consistency of the covariates
        stopifnot(all(sapply(covarList, ncol) == ncol(adjacencyMatrix))) # with the network data

        ## INITIALIZE THE SBM OBJECT ACCORDING TO THE DATA
        connectParam <- switch(model,
          "bernoulli"  = list(mean = matrix(0, 0, 0)),
          "poisson"    = list(mean = matrix(0, 0, 0)),
          "gaussian"   = list(mean = matrix(0, 0, 0), var = 1),
          "ZIgaussian" = list(mean = matrix(0, 0, 0), var = 1, p0 = 0),
        )

        super$initialize(model        = model,
                         directed     = directed,
                         nbNodes      = nrow(adjacencyMatrix),
                         blockProp    = vector("numeric", 0),
                         connectParam = connectParam,
                         dimLabels    = dimLabels,
                         covarList    = covarList)
        private$Y <- adjacencyMatrix
      },
      #--------------------------------------------
      #' @description function to perform optimization
      #' @param estimOptions a list of parameters controlling the inference algorithm and model selection. See details.
      #' @inherit estimateSimpleSBM details
      optimize = function(estimOptions = list()){

        if(private$model == 'ZIgaussian') stop("Inference not yet implemented for ZI gaussian network")

        currentOptions <- list(
          verbosity     = 3,
          plot          = TRUE,
          exploreFactor  = 1.5,
          exploreMin = 4,
          exploreMax  = Inf,
          nbBlocksRange = c(4,Inf),
          nbCores       = 2,
          fast          = TRUE
        )
        currentOptions[names(estimOptions)] <- estimOptions

        ## Transform estimOptions to a suited for blockmodels list of options
        blockmodelsOptions <- list(
          verbosity          = currentOptions$verbosity,
          plotting           = if(currentOptions$plot) character(0) else "",
          explore_min        = currentOptions$exploreMin,
          explore_max        = currentOptions$exploreMax,
          ncores             = currentOptions$nbCores,
          exploration_factor = currentOptions$exploreFactor
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
      reorder = function(){
        o <- order(private$theta$mean %*% private$pi, decreasing = TRUE)
        private$pi <- private$pi[o]
        private$theta$mean <- private$theta$mean[o, o, drop = FALSE]
        private$Z <- private$Z[, o, drop = FALSE]
      },
      #--------------------------------------------
      #' @description show method
      #' @param type character used to specify the type of SBM
      show = function(type = "Fit of a Simple Stochastic Block Model"){
        super$show(type)
        cat("* Additional fields\n")
        cat("  $probMemberships, $loglik, $ICL, $storedModels, \n")
        cat("* Additional methods \n")
        cat("  predict, fitted, $setModel, $reorder \n")
      }
    ),
    active = list(
      #' @field loglik double: approximation of the log-likelihood (variational lower bound) reached
      loglik = function(value) {private$J},
      #' @field ICL double: value of the integrated classification log-likelihood
      ICL    = function(value) {private$vICL},
      #' @field penalty double, value of the penalty term in ICL
      penalty  = function(value) {unname((self$nbConnectParam + self$nbCovariates) * log(self$nbDyads) + (self$nbBlocks-1) * log(self$nbNodes))},
      #' @field entropy double, value of the entropy due to the clustering distribution
      entropy  = function(value) {-sum(.xlogx(private$Z))},
      #' @field storedModels data.frame of all models fitted (and stored) during the optimization
      storedModels = function(value) {
        nbBlocks <- unlist(sapply(private$BMobject$memberships, function(m) ncol(m$Z)))
        nbConnectParam <- unlist(sapply(private$BMobject$model_parameters, function(param) param$n_parameters))
        U <- data.frame(
          indexModel  = 1:length(nbBlocks),
          nbParams = nbConnectParam + nbBlocks - 1,
          nbBlocks = nbBlocks,
          ICL      = private$BMobject$ICL,
          loglik   = private$BMobject$PL
          )
        U[!is.na(U$nbParams),]
      }
    )
  )
