#' Sampling of Simple SBMs
#'
#' This function samples a simple Stochastic Block Models, with various model
#' for the distribution of the edges:  Bernoulli, Poisson, or Gaussian models, and possibly with covariates
#'
#' @param nbNodes number of nodes in the network
#' @param blockProp parameters for block proportions
#' @param connectParam list of parameters for connectivity with a matrix of means 'mu' and an optional matrix of variances 'sigma2', the sizes of which must match \code{blockProp} length
#' @param model character describing the model for the relation between nodes (\code{'bernoulli'}, \code{'poisson'}, \code{'gaussian'}, ...). Default is \code{'bernoulli'}.
#' @param directed logical, directed network or not. Default is \code{FALSE}.
#' @param covariates a list of matrices with same dimension as mat describing covariates at the edge level. No covariate per Default.
#' @param covariatesParam optional vector of covariates effect. A zero length numeric vector by default.
#'
#' @return  an object with class SimpleSBM_sampler
#'
#' @examples
#' ### SIMPLE SBM - Bernoulli
#' nbNodes  <- 90
#' blockProp <- c(.5, .25, .25)      # group proportions
#' connectProb <- diag(.4, 3) + 0.05 # connectivity matrix: affiliation network
#'
#' ## Graph Sampling
#' mySampler <- sampleSimpleSBM(nbNodes, blockProp, list(mu = connectProb))
#' plot(mySampler)
#' mySampler$rBlocks()     # sample new blocks
#' mySampler$rAdjMatrix() # sample new adjacency matrix
#' plot(mySampler)
#' @export
sampleSimpleSBM <- function(nbNodes,
                            blockProp,
                            connectParam,
                            model = 'bernoulli',
                            directed = FALSE,
                            covariates = list(),
                            covariatesParam = numeric(0)) {

  mySampler <- SimpleSBM_sampler$new(model, nbNodes, directed, blockProp, connectParam, covariatesParam, covariates)
  mySampler
}
