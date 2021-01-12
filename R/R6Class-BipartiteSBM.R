#' R6 virtual class for SBM representation
#'
#' @import R6
#' @include R6Class-SBM.R
#' @export
BipartiteSBM <- # this class inherit from SBM and allow to use multipartite as a list of SimpleSBM and BipartiteSBM
  R6::R6Class(
    classname = "BipartiteSBM",
    inherit = SBM,
    private = list(
      Y   = NULL,
      tau = NULL
    ),
    public = list(
      #' @description constructor for a Simple SBM fit
      #' @param incidenceMatrix incidence (weighted) matrix
      #' @param model character (\code{'bernoulli'}, \code{'poisson'}, \code{'gaussian'})
      #' @param directed logical, directed network or not. In not, \code{incidenceMatrix} must be symmetric.
      #' @param dimLabels list of labels of each dimension (in row, in columns)
      #' @param covarList and optional list of covariates, each of whom must have the same dimension as \code{incidenceMatrix}
      initialize = function(incidenceMatrix, model,  dimLabels=list(row="row", col="col"), covarList=list()) {

        ## INITIALIZE THE Bipartite SBM OBJECT ACCORDING TO THE DATA
        super$initialize(model = model, dimension = dim(incidenceMatrix), dimLabels = dimLabels, covarList = covarList)
        private$Y <- incidenceMatrix

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
                          covarList,private$theta$p0)
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
      }
    ),
    active = list(
      #' @field memberships list of size 2: vector of memberships in row, in column.
      memberships = function(value) {lapply(private$tau, as_clustering)}
    )
  )
