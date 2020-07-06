#' Estimation of Simple SBMs
#'
#' This function performs variational inference of simple Stochastic Block Models, with various model
#' for the distribution of the edges:  Bernoulli, Poisson, or Gaussian models.
#'
#' @param netMat a matrix describing the network: either an adjacency (square) or incidence matrix with possibly weighted entries.
#' @param model character describing the model for the relation between nodes (\code{'bernoulli'}, \code{'poisson'}, \code{'gaussian'}, ...). Default is \code{'bernoulli'}.
#' @param directed logical: is the network directed or not? Only relevant when \code{type} is \code{'Simple'}. Default is \code{TRUE} if \code{netMat} is symmetric, \code{FALSE} otherwise
#' @param covariates a list of matrices with same dimension as mat describing covariates at the edge level. No covariate per Default.
#' @param estimOptions a list of parameters controlling the inference algorithm and model selection. See details.
#'
#' @details The list of parameters \code{estimOptions} essentially tunes the optimization process and the variational EM algorithm, with the following parameters
#'  \itemize{
#'  \item{"nbCores"}{integer for number of cores used. Default is 2}
#'  \item{"verbosity"}{integer for verbosity (0, 1). Default is 1}
#'  \item{"plot"}{boolean, should the ICL by dynamically plotted or not. Default is TRUE}
#'  \item{"exploreFactor"}{control the exploration of the number of groups}
#'  \item{"nbBlocksRange"}{minimal and maximal number or blocks explored}
#' }
#' @return  a list with the estimated parameters. See details...
#'
#' @examples
#' ### =======================================
#' ### SIMPLE BINARY SBM (Bernoulli model)
#'
#' ## Graph parameters & Sampling
#' nbNodes  <- 90
#' blockProp <- c(.5, .25, .25) # group proportions
#' means <- diag(.4, 3) + 0.05  # connectivity matrix: affiliation network
#' connectParam <- list(mean = means)
#' mySampler <- sampleSimpleSBM(nbNodes, blockProp, connectParam)
#' adjacencyMatrix <- mySampler$netMatrix
#'
#' ## Estimation
#' mySimpleSBM <-
#'   estimateSimpleSBM(adjacencyMatrix, 'bernoulli', estimOptions = list(plot = FALSE))
#' plot(mySimpleSBM, 'data', ordered = FALSE)
#' plot(mySimpleSBM, 'data')
#' plot(mySimpleSBM, 'expected', ordered = FALSE)
#' plot(mySimpleSBM, 'expected')
#'
#' ### =======================================
#' ### SIMPLE POISSON SBM
#'
#' ## Graph parameters & Sampling
#' nbNodes  <- 90
#' blockProp <- c(.5, .25, .25) # group proportions
#' means <- diag(15., 3) + 5    # connectivity matrix: affiliation network
#' connectParam <- list(mean = means)
#' mySampler <- sampleSimpleSBM(nbNodes, blockProp, list(mean = means), model = "poisson")
#' adjacencyMatrix <- mySampler$netMatrix
#'
#' ## Estimation
#' mySimpleSBM <- estimateSimpleSBM(adjacencyMatrix, 'poisson', estimOptions = list(plot = FALSE))
#' plot(mySimpleSBM, 'data', ordered = FALSE)
#' plot(mySimpleSBM, 'data')
#' plot(mySimpleSBM, 'expected', ordered = FALSE)
#' plot(mySimpleSBM, 'expected')
#'
#' ### =======================================
#' ### SIMPLE GAUSSIAN SBM
#'
#' ## Graph parameters & Sampling
#' nbNodes  <- 90
#' blockProp <- c(.5, .25, .25)      # group proportions
#' means <- diag(15., 3) + 5 # connectivity matrix: affiliation network
#' connectParam <- list(mean = means, var = 2)
#' mySampler <- sampleSimpleSBM(nbNodes, blockProp, connectParam, model = "gaussian")
#'
#' ## Estimation
#' mySimpleSBM <- estimateSimpleSBM(mySampler$netMatrix, 'gaussian', estimOptions = list(plot = FALSE))
#' plot(mySimpleSBM, 'data', ordered = FALSE)
#' plot(mySimpleSBM, 'data')
#' plot(mySimpleSBM, 'expected', ordered = FALSE)
#' plot(mySimpleSBM, 'expected')
#'
#' @export
estimateSimpleSBM <- function(netMat,
                              model        = 'bernoulli',
                              directed     = !isSymmetric(netMat),
                              covariates   = list(),
                              estimOptions = list()) {

  ## Set default options for estimation
  currentOptions <- list(
    verbosity     = 3,
    plot          = TRUE,
    explorFactor  = 1.5,
    nbBlocksRange = c(4,Inf),
    nbCores       = 2,
    fast          = TRUE
  )

  ## Current options are default expect for those passed by the user
  currentOptions[names(estimOptions)] <- estimOptions

  ## Construct the SBM model
  mySBM <- SimpleSBM_fit$new(netMat, model, directed, covariates)

  ## Perform optimization
  do.call(mySBM$optimize, currentOptions)

  ## reordering according to large block/large probabilities
  mySBM$reorder()

  ## Send back the SimpleSBM_fit Object
  mySBM
}

#' Estimation of Bipartite SBMs
#'
#' This function performs variational inference of bipartite Stochastic Block Models, with various model
#' for the distribution of the edges:  Bernoulli, Poisson, or Gaussian models.
#'
#' @param netMat a matrix describing the network: either an adjacency (square) or incidence matrix with possibly weighted entries.
#' @param model character describing the model for the relation between nodes (\code{'bernoulli'}, \code{'poisson'}, \code{'gaussian'}, ...). Default is \code{'bernoulli'}.
#' @param covariates a list of matrices with same dimension as mat describing covariates at the edge level. No covariate per Default.
#' @param estimOptions a list of parameters controlling the inference algorithm and model selection. See details.
#'
#' @details The list of parameters \code{estimOptions} essentially tunes the optimization process and the variational EM algorithm, with the following parameters
#'  \itemize{
#'  \item{"nbCores"}{integer for number of cores used. Default is 2}
#'  \item{"verbosity"}{integer for verbosity (0, 1). Default is 1}
#'  \item{"plot"}{boolean, should the ICL by dynamically plotted or not. Default is TRUE}
#'  \item{"exploreFactor"}{control the exploration of the number of groups}
#'  \item{"nbBlocksRange"}{minimal and maximal number or blocks explored}
#' }
#' @return  a list with the estimated parameters. See details...
#'
#' @examples
#' ### =======================================
#' ### BIPARTITE BINARY SBM (Bernoulli model)
#'
#' ## Graph parameters and Sampling
#' nbNodes <- c(60, 80)
#' blockProp <- list(c(.5, .5), c(1/3, 1/3, 1/3)) # group proportions
#' means <- matrix(runif(6), 2, 3)  # connectivity matrix
#' # In Bernoulli SBM, parameters is a list with a
#' # matrix of means 'mean' which are probabilities of connection
#' connectParam <- list(mean = means)
#' mySampler <- sampleBipartiteSBM(nbNodes, blockProp, connectParam, model = 'bernoulli')
#'
#' ## Estimation
#' myBipartiteSBM <- estimateBipartiteSBM(mySampler$netMatrix, estimOptions = list(plot = FALSE))
#' plot(myBipartiteSBM, 'expected')
#'
#' ### =======================================
#' ### BIPARTITE POISSON SBM
#'
#' ## Graph parameters & Sampling
#' nbNodes <- c(60, 80)
#' blockProp <- list(c(.5, .5), c(1/3, 1/3, 1/3)) # group proportions
#' means <- matrix(rbinom(6, 30, 0.25), 2, 3)  # connectivity matrix
#' connectParam <- list(mean = means)
#' mySampler <- sampleBipartiteSBM(nbNodes, blockProp, connectParam, model = 'poisson')
#'
#' ## Estimation
#' myBipartiteSBM <-
#'   estimateBipartiteSBM(mySampler$netMatrix, 'poisson', estimOptions = list(plot = FALSE))
#' plot(myBipartiteSBM, 'expected')
#'
#' ### =======================================
#' ### BIPARTITE GAUSSIAN SBM
#' ## Graph parameters & sampling
#' nbNodes <- c(60, 80)
#' blockProp <- list(c(.5, .5), c(1/3, 1/3, 1/3)) # group proportions
#' means <- 20 * matrix(runif(6), 2, 3)  # connectivity matrix
#' connectParam <- list(mean = means, var = 1)
#' mySampler <- sampleBipartiteSBM(nbNodes, blockProp, connectParam, model = 'gaussian')
#'
#' ## Estimation
#' myBipartiteSBM <-
#'   estimateBipartiteSBM(mySampler$netMatrix, 'gaussian', estimOptions = list(plot = FALSE))
#' plot(myBipartiteSBM, 'expected')
#'
#' @export
estimateBipartiteSBM <- function(netMat,
                                 model        = 'bernoulli',
                                 covariates   = list(),
                                 estimOptions = list()) {

  ## Set default options for estimation
  currentOptions <- list(
    verbosity     = 3,
    plot          = TRUE,
    explorFactor  = 1.5,
    nbBlocksRange = c(4,Inf),
    nbCores       = 2,
    fast          = TRUE
  )

  ## Current options are default expect for those passed by the user
  currentOptions[names(estimOptions)] <- estimOptions

  ## Construct the SBM model
  mySBM <-  BipartiteSBM_fit$new(netMat, model, covariates)

  ## Perform optimization
  do.call(mySBM$optimize, currentOptions)

  ## reordering according to large block/large probabilities
  mySBM$reorder()

  ## Send back the SimpleSBM_fit Object
  mySBM
}

