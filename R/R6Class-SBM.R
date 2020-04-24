#' R6 virtual class for SBM representation (mother class of Simple and Bipartite SBM fit and sampler)
#'
#' @import R6
SBM <- # this virtual class is the mother of all subtypes of SBM (Simple or Bipartite)
  R6::R6Class(classname = "SBM",
    ## fields for internal use (referring to the mathematical notation)
    private = list(
      model_  = NULL, # characters, the model family (bernoulli, poisson, etc.) for the edges
      dim     = NULL, # vector: number of nodes in row and in col
      pi      = NULL, # vector of parameters for block prior probabilities
      theta   = NULL, # connectivity parameters between edges
      beta    = NULL, # vector of covariates parameters
      Y       = NULL, # data matrix  (dim[1] x dim[2])
      X       = NULL  # list of covariates (list of dim[1] x dim[2] matrices)
    ),
    public = list(
      #' @description constructor for SBM
      initialize = function(model=NA, dimension=NA, blockProp=NA, connectParam=NA, covarParam=numeric(0), covarList=list()) {
        ## MODEL & PARAMETERS
        private$model_ <- model
        private$dim    <- dimension
        private$X      <- covarList
        private$pi     <- blockProp
        private$theta  <- connectParam
        private$beta   <- covarParam
      }
    ),
    ## active binding to access fields outside the class
    active = list(
      #' @field integer, the number of covariates
      nbCovariates  = function(value) {length(private$X)},
      #' @field character, the family of model for the distribution of the edges
      model         = function(value) {private$family   },
      #' @field vector of block proportions (aka prior probabilities of each block)
      blockProp     = function(value) {if (missing(value)) return(private$pi)     else private$pi     <- value},
      #' @field parameters associated to the connectivity of the SBM, e.g. matrix of inter/inter block probabilities when model is Bernoulli
      connectParam  = function(value) {if (missing(value)) return(private$theta)  else private$theta  <- value},
      #' @field vector of regression parameters assocaited with the covariates.
      covarParam    = function(value) {if (missing(value)) return(private$beta)   else private$beta   <- value},
      #' @field list of matrices of covariates
      covarList     = function(value) {if (missing(value)) return(private$X)      else private$X      <- value},
      #' @field the matrix (adjacency or incidence) encoding the network
      netMatrix     = function(value) {if (missing(value)) return(private$Y)      else private$Y      <- value}
    )
  )

#' R6 virtual class for SBM fit (mother class of Simple and Bipartite SBM fit)
#'
#' @import R6
SBM_fit <- # this virtual class is the mother of all subtypes of SBM (Simple or Bipartite)
  R6::R6Class(classname = "SBM_fit",
    inherit = SBM,
    ## fields for internal use (referring to the mathematical notation)
    private = list(
      optimizer_name = NULL,
      J              = NULL, # variational approximation of the log-likelihood
      vICL           = NULL, # variational approximation of the ICL
      tau            = NULL, # variational parameters for posterior probablility of class belonging
      Y_hat          = NULL
    ),
    public = list(
      initialize = function(data, model, covarList) {

        super$initialize(model = model, dimension = dim(data), covarList = covarList)
        private$Y <- data

        ## OPTIMIZATION FUNCTION NAME
        private$optimizer_name <- paste0("BM_", model, ifelse(self$nbCovariates > 0, "_covariates", ""))

      }
    ),
    active = list(
      #' @field size-2 vector: dimension of the network
      dimension       = function(value) {private$dim  },
      #' @field matrix -- or list of 2 matrices for Bipartite network -- of estimated probabilities for block memberships for all nodes
      probMemberships = function(value) {private$tau  },
      #' @field double: approximation of the log-likelihood (variational lower bound) reached
      loglik          = function(value) {private$J    },
      #' @field double: value of tje integrated classification loglielihood
      ICL             = function(value) {private$vICL },
      #' @field matrix of predicted value of the network
      fitted          = function(value) {private$Y_hat}
    )
  )

