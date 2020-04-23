#' @import R6
SBM <- # this virtual class is the mother of all subtypes of SBM (Simple or Bipartite)
  R6::R6Class(classname = "SBM",
    ## fields for internal use (referring to the mathematical notation)
    private = list(
      family   = NULL, # (vector of) character(s), the model family (bernoulli, poisson, etc.) for the edges
      dim      = NULL, # (list of) vector(s): number of nodes in row and in col
      pi       = NULL, # (list of) block prior probabilities
      theta    = NULL, # (list of) connectivity parameters between edges
      beta     = NULL, # vector of covariates parameters
      Y        = NULL, # (list of) data matrix(ces)  (dim[1] x dim[2] with dim[1] == dim[2] for simple SBM)
      X        = NULL  # list of covariates (list of dim[1] x dim[2] matrices)
    ),
    public = list(
      ## constructor
      initialize = function(model=NA, dimension=NA, blockProp=NA, connectParam=NA, covarParam=numeric(0), covarList=list()) {
        ## MODEL & PARAMETERS
        private$family   <- model
        private$dim      <- dimension
        private$X        <- covarList
        private$pi       <- blockProp
        private$theta    <- connectParam
        private$beta     <- covarParam
      }
    ),
    active = list(
      ## active binding to access fields outside the class
      nbCovariates  = function(value) {length(private$X)},
      model         = function(value) {private$family   },
      ## TODO --- CHECK IF WE SHOULD LET THE POSSIBILITY TO UPDATE THESE FIELDS ---
      blockProp     = function(value) {if (missing(value)) return(private$pi)     else private$pi     <- value},
      connectParam  = function(value) {if (missing(value)) return(private$theta)  else private$theta  <- value},
      covarParam    = function(value) {if (missing(value)) return(private$beta)   else private$beta   <- value},
      covarList     = function(value) {if (missing(value)) return(private$X)      else private$X      <- value},
      netData       = function(value) {if (missing(value)) return(private$Y)      else private$Y      <- value}
    )
  )

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

        if (is.list(data)) dimension <- lapply(data, dim) else dimension <- dim(data)

        super$initialize(model = model, dimension = dimension, covarList = covarList)
        private$Y <- data

        ## OPTIMIZATION FUNCTION NAME
        ### TODO -- check that all models are covered
        if (is.list(data)) {

        } else {
          private$optimizer_name <- paste0("BM_", model, ifelse(self$nbCovariates > 0, "_covariates", ""))
        }

      }
    ),
    active = list(
      dimension       = function(value) {private$dim  },
      probMemberships = function(value) {private$tau  },
      loglik          = function(value) {private$J    },
      ICL             = function(value) {private$vICL },
      fitted          = function(value) {private$Y_hat}
    )
  )
