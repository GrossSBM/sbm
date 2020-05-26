#' R6 Class definition of an Score  SBM fit
#'
#' This class is designed to give a representation and adjust an SBM fitted with blockmodels from Score matrices.
#'
#' @import R6 blockmodels
#' @include R6Class-SBM_fit.R
#' @export
ScoreSBM_fit <-
  R6::R6Class(classname = "ScoreSBM_fit",
    inherit = SBM_fit,
    private = list(
      directed_ = NULL, # is the network directed or not
      S = NULL # list of the Scores
    ),
    public = list(
      #' @description constructor for a Score SBM fit
      #' @param scores a list of square (score) matrices
      #' @param directed logical, directed network or not. In not, \code{adjacencyMatrix} must be symmetric.
      initialize = function(scores, directed) {

        ## SANITY CHECKS
        stopifnot(all.equal(nrow(scores[[1]]), ncol(scores[[1]])))  # matrix must be square
        stopifnot(isSymmetric(scores[[1]]) == !directed)                # symmetry and direction must agree

        ## INITIALIZE THE SBM OBJECT ACCORDING TO THE DATA
        adj <- matrix('NA', nrow(scores[[1]]),ncol(scores[[1]]))
        super$initialize(adj, model = 'bernoulli',covarList = list())
        private$directed_ <- directed
        private$S <- scores

      },

      #' @description function to perform optimization
      #' @param estimOptions options for the estimation
      #' @param monitoring monitoring Options
      optimize = function(estimOptions=list(),monitoring = list()) {

        resOptim <- optimizeScoreSBM(scores = private$S,
                                      directed = private$directed_,
                                      estimOptions = estimOptions,
                                      monitoring = monitoring)




    ## Exporting blockmodels output to simpleSBM_fit fields
        ind_best      <- 1
        resEstim <- resOptim[[ind_best]]
        private$J     <- resEstim$lowerBound
        private$vICL  <- resEstim$ICL
        private$adjacency_hat <- vect2Mat(resEstim$qDist$psi[,2], symmetric = !private$directed_, diag = FALSE)
        private$tau   <- resEstim$qDist$tau
        private$pi    <- resEstim$theta$blockProp
        private$theta <-  list(mu = resEstim$theta$connectParam)
        private$emiParam <- resEstim$theta$emissionParam
        names(private$emiParam$edgeParam) <- names(private$emiParam$noEdgeParam) <- c('mu', 'sigma2')

        invisible(resOptim)
      },
      #' @description prediction under the currently estimated model
      #' @return a matrix of expected values for each score and dyad
      predict = function() {
        stopifnot(is.list(covarList), self$nbCovariates == length(covarList))
        r <- private$tau %*% private$theta$mean %*% t(private$tau)
        return(r * private$emiParam$edgeParam$mu + (1 - r)*private$emiParam$noEdgeParam$mu)

      },
      #' @description show method
      #' @param type character used to specify the type of SBM
      show = function(type = "Fit of a Score  Stochastic Block Model") super$show(type)
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
      #' @field directed is the network directed or not
      directed = function(value) {private$directed_},
      #' @field scores is the list of the scores of each dyad in the network
      scores = function(value) {private$S},
      #' @field nbScores is the number of scores observed for each dyad
      nbScores = function(value) {length(private$S)}

    )
  )
