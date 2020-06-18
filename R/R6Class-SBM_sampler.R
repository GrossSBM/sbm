#' R6 class for SBM sampler
#'
#' @import R6
#' @include R6Class-SBM.R
SBM_sampler <- # Virtual call for SBM sampler (children: Simple and Bipartite SBM)
  R6::R6Class(classname = "SBM_sampler",
    inherit = SBM,
    ## fields for internal use (referring to the mathematical notation)
    private = list(
      Z = NULL, # the sampled indicator of blocks
      sampling_func = NULL # a function to sample edge values, depending on the model
    ),
    public = list(
      #' @description constructor for SBM
      #' @param model character describing the type of model
      #' @param nbNodes number of nodes in the network
      #' @param blockProp parameters for block proportions (vector of list of vectors)
      #' @param connectParam list of parameters for connectivity with a matrix of means 'mean' and an optional scalar for the variance 'var'. The dimensions of mu must match \code{blockProp} lengths
      #' @param covarParam optional vector of covariates effect
      #' @param covarList optional list of covariates data
      initialize = function(model, nbNodes, blockProp, connectParam, covarParam=numeric(0), covarList=list()) {

        super$initialize(model = model, dimension = nbNodes, blockProp = blockProp, connectParam = connectParam, covarParam = covarParam, covarList = covarList)

        ## ADDITIONAL CHECKS
        if (model == 'bernoulli') {
          stopifnot(all(connectParam$mean >= 0))
          stopifnot(all(connectParam$mean <= 1))
        }
        if (model == 'poisson') {
          stopifnot(all(connectParam$mean >= 0))
        }
        if (model == 'gaussian') {
          stopifnot(length(connectParam$var) == 1)
          stopifnot(connectParam$var > 0)
        }

        private$sampling_func <- switch(model,
            "gaussian"  = function(n, param) rnorm(n = n, mean   = param$mean, sd = sqrt(param$var)) ,
            "poisson"   = function(n, param) rpois(n = n, lambda = param$mean) ,
            "bernoulli" = function(n, param) rbinom(n = n, size = 1, prob   = param$mean),
          )
        self$rMemberships()
      },
      #' @description print/show method
      #' @param type character to tune the displayed name
      show = function(type = "Sampler for a Stochastic Block Model") {
        super$show(type)
        cat("  $indMemberships, $memberships, $expectation\n")
        cat("* S3 methods \n")
        cat("  plot, print, coef \n")
      }
    ),
    active = list(
      #' @field variance variance of each dyad under the current model
      variance = function() {if (private$model == 'gaussian') return(private$theta$var) else return(NULL) }
    )
  )


