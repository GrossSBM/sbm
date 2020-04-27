#' R6 virtual class for SBM representation (mother class of Simple and Bipartite SBM fit and sampler)
#'
#' @import R6
SBM <- # this virtual class is the mother of all subtypes of SBM (Simple or Bipartite)
  R6::R6Class(classname = "SBM",
    ## fields for internal use (referring to the mathematical notation)
    private = list(
      model   = NULL, # characters, the model name: distribution of the deges (bernoulli, poisson, gaussian) + covariates or not
      link    = NULL, # the link function like  in GLM setup
      invlink = NULL, # the inverse link function like in GLM setup
      dim     = NULL, # vector: number of nodes in row and in col
      pi      = NULL, # vector of parameters for block prior probabilities
      theta   = NULL, # connectivity parameters between edges
      beta    = NULL, # vector of covariates parameters
      Y       = NULL, # data matrix  (dim[1] x dim[2])
      X       = NULL  # list of covariates (list of dim[1] x dim[2] matrices)
    ),
    public = list(
      #' @description constructor for SBM
      #' @param model character describing the type of model
      #' @param dimension dimension of the network matrix
      #' @param blockProp parameters for block proportions (vector of list of vectors)
      #' @param connectParam matrix of parameters for connectivity
      #' @param covarParam optional vector of covariates effect
      #' @param covarList optional list of covariates data
      initialize = function(model=NA, dimension=NA, blockProp=NA, connectParam=NA, covarParam=numeric(0), covarList=list()) {

        ## MODEL & PARAMETERS
        private$model <- ifelse(length(covarList) > 0, paste0(model,"_covariates"), model)
        private$dim   <- dimension
        private$X     <- covarList
        private$pi    <- blockProp
        private$theta <- connectParam
        private$beta  <- covarParam
        private$link  <- switch(model,
                "gaussian"  = function(x) {x},
                "poisson"   = function(x) {log(x)},
                "bernoulli" = function(x) {.logit(x)},
                )
        private$invlink  <- switch(model,
                "gaussian"  = function(x) {x},
                "poisson"   = function(x) {exp(x)},
                "bernoulli" = function(x) {.logistic(x)},
                )
      }
    ),
    ## active binding to access fields outside the class
    active = list(
      #' @field dimension size-2 vector: dimension of the network
      dimension       = function(value) {private$dim  },
      #' @field modelName character, the family of model for the distribution of the edges
      modelName    = function(value) {length(private$model)},
      #' @field nbCovariates integer, the number of covariates
      nbCovariates = function(value) {length(private$X)},
      #' @field blockProp vector of block proportions (aka prior probabilities of each block)
      blockProp     = function(value) {if (missing(value)) return(private$pi)     else private$pi     <- value},
      #' @field connectParam parameters associated to the connectivity of the SBM, e.g. matrix of inter/inter block probabilities when model is Bernoulli
      connectParam = function(value) {if (missing(value)) return(private$theta)  else private$theta  <- value},
      #' @field covarParam vector of regression parameters assocaited with the covariates.
      covarParam   = function(value) {if (missing(value)) return(private$beta)   else private$beta   <- value},
      #' @field covarList list of matrices of covariates
      covarList    = function(value) {if (missing(value)) return(private$X)      else private$X      <- value},
      #' @field covarEffect effect of covariates
      covarEffect = function(value) {roundProduct(simplify2array(private$X), private$beta)},
      #' @field netMatrix the matrix (adjacency or incidence) encoding the network
      netMatrix    = function(value) {if (missing(value)) return(private$Y)      else private$Y      <- value}
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
      J              = NULL, # variational approximation of the log-likelihood
      vICL           = NULL, # variational approximation of the ICL
      tau            = NULL, # variational parameters for posterior probablility of class belonging
      Y_hat          = NULL
    ),
    public = list(
      #' @description constructor for SBM fit
      #' @param data the data matrix of the networj
      #' @param model character describing the type of model
      #' @param covarList optional list of matrices for covariates
      initialize = function(data, model, covarList) {
        super$initialize(model = model, dimension = dim(data), covarList = covarList)
        private$Y <- data
      }
    ),
    active = list(
      #' @field probMemberships matrix -- or list of 2 matrices for Bipartite network -- of estimated probabilities for block memberships for all nodes
      probMemberships = function(value) {private$tau  },
      #' @field loglik double: approximation of the log-likelihood (variational lower bound) reached
      loglik          = function(value) {private$J    },
      #' @field ICL double: value of tje integrated classification loglielihood
      ICL             = function(value) {private$vICL },
      #' @field fitted matrix of predicted value of the network
      fitted          = function(value) {private$Y_hat}
    )
  )

## PUBLIC S3 METHODS FOR SBM
## =========================================================================================

## Auxiliary functions to check the given class of an objet
is_SBM <- function(Robject) {inherits(Robject, "SBM")}
is_SBM_fit <- function(Robject) {inherits(Robject, "SBM_fit")}

#' Extract model coefficients
#'
#' Extracts model coefficients from objects with class \code{\link[=SBM]{SBM}} and children (\code{\link[=SimpleSBM]{SimpleSBM}},
#' \code{\link[=SimpleSBM]{BipartiteSBM}})
#'
#' @param object an R6 object inheriting from class SBM_fit (like SimpleSBM_fit or BipartiteSBM_fit)
#' @param type type of parameter that should be extracted. Either 'memberships' for \deqn{\pi}, 'connectivity' for \deqn{\theta},
#'  or "covariates" for \deqn{\beta}. Default is 'connectivity'.
#' @param ... additional parameters for S3 compatibility. Not used
#' @return vector or list of parameters.
#' @export
coef.SBM <- function(object, type = c( 'connectivity','membership', 'covariates'), ...) {
  stopifnot(is_SBM(object))
  switch(match.arg(type),
         membership   = object$blockParam,
         connectivity = object$connectParam,
         covariates   = object$covarParam)
}

#' @importFrom stats fitted
#' @export
predict.SBM_fit <- function(object, ...) {
  stopifnot(is_SBM(object))
  object$connectProb
}

#' @export
summary.SBM <- function(object, ...) {
  stopifnot(is_SBM(object))
  object$show()
}

