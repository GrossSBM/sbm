#' Estimation of Simple SBMs
#'
#' This function performs variational inference of simple Stochastic Block Models, with various model
#' for the distribution of the edges:  Bernoulli, Poisson, or Gaussian models.
#'
#' @param netMat a matrix describing the network: either an adjacency (square) or incidence matrix with possibly weighted entries.
#' @param model character describing the model for the relation between nodes (\code{'bernoulli'}, \code{'poisson'}, \code{'gaussian'}, ...). Default is \code{'bernoulli'}.
#' @param directed logical: is the network directed or not? Only relevant when \code{type} is \code{'Simple'}. Default is \code{TRUE} if \code{netMat} is symmetric, \code{FALSE} otherwise
#' @param dimLabels an optional list of labels for each dimension (in row, in column)
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
                              dimLabels    = list(row = "rowLabel", col = "colLabel"),
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
  mySBM <- SimpleSBM_fit$new(netMat, model, directed, dimLabels, covariates)

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
#' @param dimLabels an optional list of labels for each dimension (in row, in column)
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
                                 dimLabels    = list(row = "rowLabel", col = "colLabel"),
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
  mySBM <-  BipartiteSBM_fit$new(netMat, model, dimLabels, covariates)

  ## Perform optimization
  do.call(mySBM$optimize, currentOptions)

  ## reordering according to large block/large probabilities
  mySBM$reorder()

  ## Send back the SimpleSBM_fit Object
  mySBM
}


#' Estimation for multipartite SBM
#'
#' @param ldefinedNet list of networks that were defined by the \code{defineSBM} function
#' @param estimOptions options for the inference procedure
#'
#' @return a MultipartiteSBM_fit object with the estimated parameters and the blocks in each Functional Group
#' @export
#'
#' @examples
#' #' # About the Functional Groups (FG)
#' nbFunctionalGroups <- 3  #number of functional groups
#' nbBlocks  <- c(3,2,2) #number of clusters in each functional group
#' nbNodes <-  c(100,50,40)
#' blockProp <- vector("list", 3)  # parameters of clustering in each functional group
#' blockProp[[1]] <- c(0.4,0.3,0.3) # in Functional Group 1
#' blockProp[[2]] <- c(0.6,0.4) # in Functional Group 2
#' blockProp[[3]]  <- c(0.6,0.4) # in Functional Group 3
#' # About the interactions between the FG
#' archiMultipartite  <-  rbind(c(1,2),c(2,3),c(2,2),c(1,3)) # architecture of the various networks (FG interaction : 1 with  2, 2 wih 3, 1 with 3 and interactions inside FG 2. )
#' model <- c('bernoulli','poisson','bernoulli','gaussian') # type of distribution in each network
#' directed <- c( NA, NA  ,  FALSE , NA) # for each network : directed or not (not required for an interaction wetween two different FG)
#' connectParam <- list()
#' E <- archiMultipartite
#' connectParam[[1]] <- list(mean = matrix(rbeta(nbBlocks[E[1,1]] * nbBlocks[E[1,2]],1,1 ),nrow = nbBlocks[E[1,1]], ncol = nbBlocks[E[1,2]] ))
#' connectParam[[2]] <- list(mean  =  matrix(rgamma(nbBlocks[E[2,1]] * nbBlocks[E[2,2]],7.5,0.01 ),nrow = nbBlocks[E[2,1]], ncol = nbBlocks[E[2,2]]))
#' connectParam[[3]] <- list(mean  =  matrix(rbeta(nbBlocks[E[3,1]] * nbBlocks[E[3,2]],0.9,0.0 ), nrow = nbBlocks[E[3,1]], ncol = nbBlocks[E[3,2]]))
#' connectParam[[3]]$mean <-  0.5*(connectParam[[3]]$mean + t(connectParam[[3]]$mean)) # symetrisation for network 3
#' connectParam[[4]] <- list(mean = matrix(rnorm(nbBlocks[E[4,1]] * nbBlocks[E[4,2]],7.5,10 ), nrow = nbBlocks[E[4,1]], ncol = nbBlocks[E[4,2]]))
#' connectParam[[4]]$var <- matrix(rgamma(nbBlocks[E[4,1]] * nbBlocks[E[4,2]],7.5,0.1 ), nrow = nbBlocks[E[4,1]], ncol = nbBlocks[E[4,2]])
#' ## Graph Sampling
#' mySampleMSBM <- sampleMultipartiteSBM(nbNodes, blockProp, archiMultipartite, connectParam, model, directed, dimLabels = as.list(c('A','B','C')))
#' listSBM <- mySampleMSBM$listSBM
#' estimOptions = list(initBM = FALSE)
#' myMSBM <- estimateMultipartiteSBM(listSBM,estimOptions)
estimateMultipartiteSBM <- function(listSBM,
                                    estimOptions = list())
{

  myMSBM <- MultipartiteSBM_fit$new(listSBM)


  currentOptions <- list(
    verbosity     = 3,
    nbBlocksRange = lapply(1:myMSBM$nbLabels,function(l){c(1,10)}),
    nbCores       = 2,
    maxiterVE     = NULL,
    maxiterVEM    = NULL,
    initBM = TRUE
    )
  names(currentOptions$nbBlocksRange) <- myMSBM$dimLabels
  ## Current options are default expect for those passed by the user
  currentOptions[names(estimOptions)] <- estimOptions

  myMSBM$optimize(currentOptions)

  myMSBM
}
