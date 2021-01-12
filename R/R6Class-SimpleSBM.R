#' R6 virtual class for SBM representation (mother class of Simple and Bipartite SBM fit and sampler)
#'
#' @import R6
#' @include R6Class-SBM.R
#' @export
SimpleSBM <- # this class inherit from SBM and allow to use multipartite as a list of SimpleSBM and BipartiteSBM
  R6::R6Class(
    classname = "SimpleSBM",
    inherit = SBM,
    private = list(
      directed_ = NULL,
      Y         = NULL,
      tau       = NULL
    ),
    public = list(
      #' @description constructor for a Simple SBM fit
      #' @param adjacencyMatrix square (weighted) matrix
      #' @param model character (\code{'bernoulli'}, \code{'poisson'}, \code{'gaussian'})
      #' @param directed logical, directed network or not. In not, \code{adjacencyMatrix} must be symmetric.
      #' @param dimLabels character labels of the nodes (in row, in columns)
      #' @param covarList and optional list of covariates, each of whom must have the same dimension as \code{adjacencyMatrix}
      initialize = function(adjacencyMatrix, model, directed, dimLabels= "node", covarList=list()) {

        ## SANITY CHECKS
        stopifnot(all.equal(nrow(adjacencyMatrix), ncol(adjacencyMatrix)))  # matrix must be square
        stopifnot(isSymmetric(adjacencyMatrix) == !directed)                # symmetry and direction must agree

        ## INITIALIZE THE SBM OBJECT ACCORDING TO THE DATA
        super$initialize(model = model, dimension = dim(adjacencyMatrix), dimLabels = dimLabels, covarList = covarList)
        private$Y <- adjacencyMatrix
        private$directed_ <- directed
      },
      #' @description prediction under the currently parameters
      #' @param covarList a list of covariates. By default, we use the covariates with which the model was estimated
      #' @return a matrix of expected values for each dyad
      predict = function(covarList = self$covarList) {
        mu <- predict_sbm(self$nbNodes,
                          self$nbCovariates,
                          private$link,
                          private$invlink,
                          private$tau,
                          private$theta$mean,
                          self$covarEffect,
                          covarList)
        mu
      },
      #' @description permute group labels by order of decreasing probability
      reorder = function(){
        o <- order_sbm(private$theta$mean,private$pi)
        private$pi <- private$pi[o]
        private$theta$mean <- private$theta$mean[o,o]
        private$tau <- private$tau[, o, drop = FALSE]
      }

    ),
    active = list(
      #' @field memberships vector of clustering
      memberships = function(value) {as_clustering(private$tau)},
      #' @field directed is the network directed or not
      directed = function(value) {private$directed_}
    )
  )

