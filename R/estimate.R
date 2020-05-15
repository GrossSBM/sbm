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
#'  \item{"nbCores"}{integer for number of cores used. Default is 1.}
#'  \item{"verbosity"}{integer for verbosity (0, 1). Default is 1.}
#'  \item{"exploreFactor"}{control the exploration of the number of groups}
#'  \item{"nbBlocksRange"}{minimal and maximal number or blocks explored}
#' }
#' @return  a list with the estimated parameters. See details...
#'
#' @importFrom corrplot corrplot
#'
#' @examples
#' ### SIMPLE SBM
#'
#' ## Graph parameters
#' nbNodes  <- 90
#' blockProp <- c(.5, .25, .25)      # group proportions
#' connectProb <- diag(.4, 3) + 0.05 # connectivity matrix: affiliation network
#'
#' ## Sampling
#' mySampler <- sampleSimpleSBM(nbNodes, blockProp, list(mu = connectProb))
#' adjacencyMatrix <- mySampler$netMatrix
#'
#' ## Estimation
#' mySimpleSBM <- estimateSimpleSBM(adjacencyMatrix)
#' par(mfrow = c(2,2))
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
    nbCores       = 1
  )

  ## Current options are default expect for those passed by the user
  currentOptions[names(estimOptions)] <- estimOptions

  ## Construct the SBM model
  mySBM <- SimpleSBM_fit$new(netMat, model, directed, covariates)

  ## Perform optimization
  do.call(mySBM$optimize, currentOptions)

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
#'  \item{"nbCores"}{integer for number of cores used. Default is 1.}
#'  \item{"verbosity"}{integer for verbosity (0, 1). Default is 1.}
#'  \item{"exploreFactor"}{control the exploration of the number of groups}
#'  \item{"nbBlocksRange"}{minimal and maximal number or blocks explored}
#' }
#' @return  a list with the estimated parameters. See details...
#'
#' @importFrom corrplot corrplot
#'
#' @examples
#' ### BIPARTITE SBM
#'
#' ## Graph parameters
#' nbNodes <- c(100, 120)
#' blockProp <- list(c(.5, .5), c(1/3,1/3,1/3)) # group proportions
#' connectProb <- matrix(runif(6), 2, 3)
#'
#' ## Graph Sampling
#' mySampler <- sampleBipartiteSBM(nbNodes, blockProp, list(mu = connectProb))
#'
#' ## Estimation
#' myBipartiteSBM <- estimateBipartiteSBM(mySampler$netMatrix)
#' plot(myBipartiteSBM, 'expected')
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
    nbCores       = 1
  )

  ## Current options are default expect for those passed by the user
  currentOptions[names(estimOptions)] <- estimOptions

  ## Construct the SBM model
  mySBM <-  BipartiteSBM_fit$new(netMat, model, covariates)

  ## Perform optimization
  do.call(mySBM$optimize, currentOptions)

  ## Send back the SimpleSBM_fit Object
  mySBM
}

