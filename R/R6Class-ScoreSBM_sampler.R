#' R6 class for Score SBM sampler
#'
#' @import R6
#' @include R6Class-SimpleSBM_sampler.R
#' @export
ScoreSBM_sampler <-
  R6::R6Class(classname = "ScoreSBM_sampler",
   inherit = SimpleSBM_sampler,
   ## fields for internal use (referring to the mathematical notation)
    private = list(
      d = NULL,   # number of Scores
      emiParam = NULL, #  parameters of the Score distribution
      S  = NULL
    ),
    public = list(
      #' @description constructor for ScoreSBM
      #' @param nbNodes number of nodes in the network
      #' @param directed logical, directed network or not.
      #' @param blockProp parameters for block proportions (vector of list of vectors)
      #' @param connectParam list of parameters for connectivity with a matrix of means 'mu' and an optional scalar for the variance 'sigma2'. The size of mu must match \code{blockProp} length
      #' @param emissionParam parameters of the emission of the Noisy SBM. List of two terms : noEdgeParam and edgeParam. Each element is a list of mean (vector or size nbScores) and variance matrix (var: square matrix of size nbScores)
      #' @param seed integer to set seed
      initialize = function(nbNodes, directed, blockProp, connectParam, emissionParam, seed = NULL) {

        set.seed(seed)
        ## ADDITIONAL SANITY CHECKS
        dims <- c(length(emissionParam$edgeParam$mu),dim(emissionParam$noEdgeParam$sigma2),dim(emissionParam$edgeParam$sigma2))
        stopifnot(all(length(emissionParam$noEdgeParam$mu) == dims))
        stopifnot(isSymmetric(emissionParam$edgeParam$sigma2))
        stopifnot(isSymmetric(emissionParam$noEdgeParam$sigma2))
        stopifnot(min(eigen(emissionParam$edgeParam$sigma2,symmetric = TRUE)$value) > 0)
        stopifnot(min(eigen(emissionParam$noEdgeParam$sigma2,symmetric = TRUE)$value) > 0)

        super$initialize('bernoulli', nbNodes, directed, blockProp, connectParam)
        private$d <- length(emissionParam$noEdgeParam$mu)
        private$emiParam <- emissionParam
        private$S <- NULL
        self$rScores()
      },
      #' @description a method to sample the score matrices
      #' @return nothing (sampled Scores are stored in the current object)
      rScores = function() {
          YArray <- replicate(private$d, private$Y, simplify = "array")
          S1 <- mvtnorm::rmvnorm(self$nbNodes**2,private$emiParam$edgeParam$mu,private$emiParam$edgeParam$sigma2)
          S0 <- mvtnorm::rmvnorm(self$nbNodes**2,private$emiParam$noEdgeParam$mu,private$emiParam$noEdgeParam$sigma2)
          SArray <- array(S1, c(self$nbNodes, self$nbNodes, private$d)) * YArray + array(S0, c(self$nbNodes, self$nbNodes, private$d)) * (1 - YArray)
          S <- lapply(1:private$d, function(i) {U <- SArray[, , i]; diag(U) <- NA; return(U)})
          if (!private$directed_) {S <- lapply(1:private$d, function(i) {0.5 * (S[[i]] + t(S[[i]]))})}
          private$S <- S
      },
      #' @description show method
      #' @param type character used to specify the type of SBM
      show = function(type = "Sampler for a Score Stochastic Block Model") {
        cat(type, "--\n")

        cat("=====================================================================\n")

        cat("Dimension = (", self$dimension, ") - (",
            self$nbBlocks, ") blocks and", "(" ,private$d, ")",  "score(s).\n")
        cat("=====================================================================\n")

        cat("* Useful fields \n")
        cat("  $nbNodes, $nbBlocks, $nbDyads\n")
        cat("  $blockProp, $connectParam, $nbScores, $emissionParam \n")
        cat("  $memberships, $netMatrix, $scores, \n")

        cat("* R6 methods \n")
        cat("  $rMemberships(), $rAdjacency(),$rScores() \n")

       }
    ),
    active = list(
    #' @field nbScores number of scores
    nbScores = function(value) {private$d},
    #' @field emissionParam parameters of the emission distribution of the scores
    emissionParam = function(value) {private$emiParam},
    #' @field scores a list of length nbScores, each matrix correspondonding to a score for each dyad
    scores = function(value) {private$S}
    )
  )
