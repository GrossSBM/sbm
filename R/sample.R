#' Sampling of Simple SBMs
#'
#' This function samples a simple Stochastic Block Models, with various model
#' for the distribution of the edges:  Bernoulli, Poisson, or Gaussian models, and possibly with covariates
#'
#' @param nbNodes number of nodes in the network
#' @param blockProp parameters for block proportions
#' @param connectParam list of parameters for connectivity with a matrix of means 'mean' and an optional matrix of variances 'var', the sizes of which must match \code{blockProp} length
#' @param model character describing the model for the relation between nodes (\code{'bernoulli'}, \code{'poisson'}, \code{'gaussian'}, ...). Default is \code{'bernoulli'}.
#' @param directed logical, directed network or not. Default is \code{FALSE}.
#' @param covariates a list of matrices with same dimension as mat describing covariates at the edge level. No covariate per Default.
#' @param covariatesParam optional vector of covariates effect. A zero length numeric vector by default.
#'
#' @return  an object with class \code{\link{SimpleSBM_sampler}}
#'
#' @examples
#' ### =======================================
#' ### SIMPLE BINARY SBM (Bernoulli model)
#' ## Graph parameters
#' nbNodes  <- 90
#' blockProp <- c(.5, .25, .25) # group proportions
#' means <- diag(.4, 3) + 0.05  # connectivity matrix: affiliation network
#' # In Bernoulli SBM, parameters is a list with a
#' # matrix of means 'mean' which are probabilities of connection
#' connectParam <- list(mean = means)
#'
#' ## Graph Sampling
#' mySampler <- sampleSimpleSBM(nbNodes, blockProp, connectParam, model = 'bernoulli')
#' plot(mySampler)
#' mySampler$rMemberships() # sample new memberships
#' mySampler$rAdjacency()   # sample new adjacency matrix
#' plot(mySampler)
#' hist(mySampler$netMatrix)
#'
#' ### =======================================
#' ### SIMPLE POISSON SBM
#' ## Graph parameters
#' nbNodes  <- 90
#' blockProp <- c(.5, .25, .25) # group proportions
#' means <- diag(15., 3) + 5    # connectivity matrix: affiliation network
#' # In Poisson SBM, parameters is a list with
#' # a matrix of means 'mean' which are a mean integer value taken by edges
#' connectParam <- list(mean = means)
#'
#' ## Graph Sampling
#' mySampler <- sampleSimpleSBM(nbNodes, blockProp, list(mean = means), model = "poisson")
#' plot(mySampler)
#' hist(mySampler$netMatrix)
#'
#' ### =======================================
#' ### SIMPLE GAUSSIAN SBM
#' ## Graph parameters
#' nbNodes  <- 90
#' blockProp <- c(.5, .25, .25)      # group proportions
#' means <- diag(15., 3) + 5 # connectivity matrix: affiliation network
#' # In Gaussian SBM, parameters is a list with
#' # a matrix of means 'mean' and a matrix of variances 'var'
#' connectParam <- list(mean = means, var = 2)
#'
#' ## Graph Sampling
#' mySampler <- sampleSimpleSBM(nbNodes, blockProp, connectParam, model = "gaussian")
#' plot(mySampler)
#' hist(mySampler$netMatrix)
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

#' Sampling of Bipartite SBMs
#'
#' This function samples a simple Stochastic Block Models, with various model
#' for the distribution of the edges:  Bernoulli, Poisson, or Gaussian models, and possibly with covariates
#'
#' @param nbNodes number of nodes in the network
#' @param blockProp parameters for block proportions: list of size two with row and column block proportions
#' @param connectParam list of parameters for connectivity with a matrix of means 'mean' and an optional matrix of variances 'var', the sizes of which must match \code{blockProp} length (in row, respectively in column)
#' @param model character describing the model for the relation between nodes (\code{'bernoulli'}, \code{'poisson'}, \code{'gaussian'}, ...). Default is \code{'bernoulli'}.
#' @param covariates a list of matrices with same dimension as mat describing covariates at the edge level. No covariate per Default.
#' @param covariatesParam optional vector of covariates effect. A zero length numeric vector by default.
#'
#' @return an object with class \code{\link{BipartiteSBM_sampler}}
#'
#' @examples
#' ### =======================================
#' ### BIPARTITE BERNOULLI SBM
#' ## Graph parameters
#' nbNodes <- c(100, 120)
#' blockProp <- list(c(.5, .5), c(1/3, 1/3, 1/3)) # group proportions
#' means <- matrix(runif(6), 2, 3)  # connectivity matrix
#' # In Bernoulli SBM, parameters is a list with
#' # a matrix of means 'mean' which are probabilities of connection
#' connectParam <- list(mean = means)
#'
#' ## Graph Sampling
#' mySampler <- sampleBipartiteSBM(nbNodes, blockProp, connectParam, model = 'bernoulli')
#' plot(mySampler)
#' mySampler$rMemberships() # sample new memberships
#' mySampler$rIncidence()   # sample new incidence matrix
#' plot(mySampler)
#' hist(mySampler$netMatrix)
#'
#' ### =======================================
#' ### BIPARTITE POISSON SBM
#' ## Graph parameters
#' nbNodes <- c(100, 120)
#' blockProp <- list(c(.5, .5), c(1/3, 1/3, 1/3)) # group proportions
#' means <- matrix(rbinom(6, 30, 0.25), 2, 3)  # connectivity matrix
#' # In Poisson SBM, parameters is a list with a matrix of
#' # means 'mean' which are a mean integer value taken by edges
#' connectParam <- list(mean = means)
#'
#' ## Graph Sampling
#' mySampler <- sampleBipartiteSBM(nbNodes, blockProp, connectParam, model = 'poisson')
#' plot(mySampler)
#' hist(mySampler$netMatrix)
#'
#' ### =======================================
#' ### BIPARTITE GAUSSIAN SBM
#' ## Graph parameters
#' nbNodes <- c(100, 120)
#' blockProp <- list(c(.5, .5), c(1/3, 1/3, 1/3)) # group proportions
#' means <- 20 * matrix(runif(6), 2, 3)  # connectivity matrix
#' # In Gaussian SBM, parameters is a list with a matrix
#' # of means 'mean' and a matrix of variances 'var'
#' connectParam <- list(mean = means, var = 1)
#'
#' ## Graph Sampling
#' mySampler <- sampleBipartiteSBM(nbNodes, blockProp, connectParam, model = 'gaussian')
#' plot(mySampler)
#' hist(mySampler$netMatrix)
#'
#' @export
sampleBipartiteSBM <- function(nbNodes,
                            blockProp,
                            connectParam,
                            model = 'bernoulli',
                            covariates = list(),
                            covariatesParam = numeric(0)) {

  mySampler <- BipartiteSBM_sampler$new(model, nbNodes, blockProp, connectParam, covariatesParam, covariates)
  mySampler
}
