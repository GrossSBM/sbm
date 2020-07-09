#' R6 virtual class for SBM fit (mother class of Simple and Bipartite SBM fit)
#'
#' @import R6
#' @include R6Class-SBM.R
SBM_fit <- # this virtual class is the mother of all subtypes of SBM (Simple or Bipartite)
  R6::R6Class(classname = "SBM_fit",
    inherit = SBM,
    ## fields for internal use (referring to the mathematical notation)
    private = list(
      J              = NULL, # variational approximation of the log-likelihood
      vICL           = NULL, # variational approximation of the ICL
      tau            = NULL, # variational parameters for posterior probability of class belonging
      BMobject       = NULL,
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
          "gaussian_covariates"       = list(mean = parameters$mu, var = parameters$sigma2)
        )
      }
    ),
    public = list(
      #' @description constructor for SBM fit
      #' @param data the data matrix of the network
      #' @param model character describing the type of model
      #' @param covarList optional list of matrices for covariates
      initialize = function(data, model, covarList) {
        super$initialize(model = model, dimension = dim(data), covarList = covarList)
        private$Y <- data
      },
      #' @description print/show method
      #' @param type character to tune the displayed name
      show = function(type = "Fit of a Stochastic Block Model") {
        super$show(type)
        cat("  $probMemberships, $memberships, $loglik, $ICL, $storedModels, $setModel \n")
        cat("* S3 methods \n")
        cat("  plot, print, coef, predict, fitted \n")
      },
      #' @description method to select a specific model among the ones fitted during the optimization.
      #'  Fields of the current SBM_fit will be updated accordingly.
      #' @param index integer, the index of the model to be selected (row number in storedModels)
      setModel = function(index) {
        stopifnot(!is.null(private$BMobject))
        stopifnot(index %in% seq.int(nrow(self$storedModels)))
        private$import_from_BM(index)
        self$reorder()
      }
    ),
    active = list(
      #' @field probMemberships matrix -- or list of 2 matrices for Bipartite network -- of estimated probabilities for block memberships for all nodes
      probMemberships = function(value) {private$tau  },
      #' @field loglik double: approximation of the log-likelihood (variational lower bound) reached
      loglik          = function(value) {private$J    },
      #' @field ICL double: value of the integrated classification log-likelihood
      ICL             = function(value) {private$vICL },
      #' @field expectation expected values of connection under the currently adjusted model
      expectation = function() {self$predict()},
      #' @field fitted matrix of predicted value of the network
      fitted          = function(value) {self$predict()}
    )
  )

# ========================================================================================
# PUBLIC S3 METHODS FOR SBM_fit

## Auxiliary function to check the given class of an objet
is_SBM_fit <- function(Robject) {inherits(Robject, "SBM_fit")}

#' Extract model fitted values
#' Extracts fitted values for object with class \code{\link[=SBM_fit]{SBM_fit}} and children (\code{\link[=SimpleSBM_fit]{SimpleSBM_fit}},
#' \code{\link[=BipartiteSBM_fit]{BipartiteSBM_fit}})
#'
#' @param object an R6 object inheriting from class SBM_fit (like SimpleSBM_fit or BipartiteSBM_fit)
#' @param ... additional parameters for S3 compatibility. Not used
#' @return a matrix of expected fitted values for each dyad
#' @importFrom stats fitted
#' @export
fitted.SBM_fit <- function(object,  ...) {
  stopifnot(is_SBM_fit(object))
  object$fitted
}

