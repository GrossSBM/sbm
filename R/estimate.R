#' Estimation of Simple SBMs
#'
#' This function performs variational inference of simple Stochastic Block Models, with various model
#' for the distribution of the edges:  Bernoulli, Poisson, or Gaussian models.
#'
#' @param netMat a matrix describing the network: either an adjacency (square) or incidence matrix with possibly weighted entries.
#' @param model character describing the model for the relation between nodes (\code{'bernoulli'}, \code{'poisson'}, \code{'gaussian'}, ...). Default is \code{'bernoulli'}.
#' @param directed logical: is the network directed or not? Only relevant when \code{type} is \code{'Simple'}. Default is \code{TRUE} if \code{netMat} is symmetric, \code{FALSE} otherwise
#' @param dimLabels an optional label for referring to the nodes
#' @param covariates a list of matrices with same dimension as mat describing covariates at the edge level. No covariate per Default.
#' @param estimOptions a list of parameters controlling the inference algorithm and model selection. See details.
#'
#' @details The list of parameters \code{estimOptions} essentially tunes the optimization process and the variational EM algorithm, with the following parameters
#' * "nbCores integer for number of cores used. Default is 2
#' * "verbosity" integer for verbosity (0, 1). Default is 1
#' * "plot" boolean, should the ICL by dynamically plotted or not. Default is TRUE
#' * "exploreFactor" control the exploration of the number of groups
#' * "exploreMin" explore at least until exploreMin even if the exploration factor rule is achieved. Default 4. See the package blockmodels for details.
#' * "exploreMax" Stop exploration at exploreMax  even if the exploration factor rule is not achieved. Default Inf. See the package blockmodels for details.
#' * "nbBlocksRange" minimal and maximal number or blocks explored
#' * "fast" logical: should approximation be used for Bernoulli model with covariates. Default to \code{TRUE}
#'
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
#' adjacencyMatrix <- mySampler$networkData
#'
#' ## Estimation
#' mySimpleSBM <-
#'   estimateSimpleSBM(adjacencyMatrix, 'bernoulli', estimOptions = list(plot = FALSE))
#' plot(mySimpleSBM, 'data', ordered = FALSE)
#' plot(mySimpleSBM, 'data')
#' plot(mySimpleSBM, 'expected', ordered = FALSE)
#' plot(mySimpleSBM, 'expected')
#' plot(mySimpleSBM, 'meso')
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
#' adjacencyMatrix <- mySampler$networkData
#'
#' ## Estimation
#' mySimpleSBM <- estimateSimpleSBM(adjacencyMatrix, 'poisson',
#'    estimOptions = list(plot = FALSE))
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
#' mySimpleSBM <-
#'    estimateSimpleSBM(mySampler$networkData, 'gaussian', estimOptions = list(plot = FALSE))
#' plot(mySimpleSBM, 'data', ordered = FALSE)
#' plot(mySimpleSBM, 'data')
#' plot(mySimpleSBM, 'expected', ordered = FALSE)
#' plot(mySimpleSBM, 'expected')
#'
#' @export
estimateSimpleSBM <- function(netMat,
                              model        = 'bernoulli',
                              directed     = !isSymmetric(netMat),
                              dimLabels    = c("node"),
                              covariates   = list(),
                              estimOptions = list()) {


  ## Set default options for estimation
  if (!is.null(estimOptions$verbosity)){
    if ((estimOptions$verbosity == 0)&(is.null(estimOptions$plot))){estimOptions$plot   = FALSE}
  }

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
  mySBM$optimize(currentOptions)

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
#' @param dimLabels an optional vector of labels for each dimension (in row, in column)
#' @param covariates a list of matrices with same dimension as mat describing covariates at the edge level. No covariate per Default.
#' @param estimOptions a list of parameters controlling the inference algorithm and model selection. See details.
#'
#' @inherit estimateSimpleSBM details
#'
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
#' myBipartiteSBM <- estimateBipartiteSBM(mySampler$networkData, estimOptions = list(plot = FALSE))
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
#'   estimateBipartiteSBM(mySampler$networkData, 'poisson', estimOptions = list(plot = FALSE))
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
#'   estimateBipartiteSBM(mySampler$networkData, 'gaussian', estimOptions = list(plot = FALSE))
#' plot(myBipartiteSBM, 'expected')
#'
#' @export
estimateBipartiteSBM <- function(netMat,
                                 model        = 'bernoulli',
                                 dimLabels    = c(row = "row", col = "col"),
                                 covariates   = list(),
                                 estimOptions = list()) {

  ## Set default options for estimation
  if (!is.null(estimOptions$verbosity)){
    if ((estimOptions$verbosity == 0)&(is.null(estimOptions$plot))){estimOptions$plot   = FALSE}
  }
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
  mySBM$optimize(currentOptions)

  ## reordering according to large block/large probabilities
  mySBM$reorder()

  ## Send back the SimpleSBM_fit Object
  mySBM
}

#-----------------------------------------------------------------
#' Estimation for multipartite SBM
#'
#' @param listSBM list of networks that were defined by the \code{defineSBM} function
#' @param estimOptions options for the inference procedure
#' @details The list of parameters \code{estimOptions} essentially tunes the optimization process and the variational EM algorithm, with the following parameters
#' * "nbCores" integer for number of cores used.  Default is 2
#' * "verbosity" integer for verbosity (0, 1). Default is 1
#' * "nbBlocksRange" List of length the number of functional groups, each element supplying the minimal and maximal number of blocks to be explored. The names of the list must be the names of the functional groups.  Default value is from 1 to 10)
#' * "initBM" Boolean. True if using simple and bipartite SBM as initialisations. Default value  = TRUE
#' * "maxiterVEM" Number of max. number of iterations in  the VEM. Default value  = 100
#' * "maxiterVE" Number of max. number of iterations in  the VE. Default value  = 100
#'
#' @return a MultipartiteSBM_fit object with the estimated parameters and the blocks in each Functional Group
#' @export
#'
#' @examples
#' \dontrun{
#' ## About the Parts/Functional Groups (FG)
#' blockProp <- list(c(0.16 ,0.40 ,0.44),c(0.3,0.7)) # prop of blocks in each FG
#' archiMultipartite <-  rbind(c(1,2),c(2,2),c(1,1)) # architecture of the multipartite net.
#' nbNodes <- c(60,50)
#' ## About the connection matrices
#' directed <- c(NA, TRUE, FALSE) # type of each network
#' model <- c('gaussian','bernoulli','poisson')
#' C1 <-
#'  list(mean = matrix(c(6.1, 8.9, 6.6, 9.8, 2.6, 1.0), 3, 2),
#'       var  = matrix(c(1.6, 1.6, 1.8, 1.7 ,2.3, 1.5),3, 2))
#' C2 <- list(mean = matrix(c(0.7,1.0, 0.4, 0.6),2, 2))
#' m3 <- matrix(c(2.5, 2.6 ,2.2 ,2.2, 2.7 ,3.0 ,3.6, 3.5, 3.3),3,3 )
#' C3 <- list(mean = .5 * (m3 + t(m3)))
#' connectParam <- list(C1, C2, C3)
#' ## Graph Sampling
#' mySampleMSBM <- sampleMultipartiteSBM(nbNodes, blockProp,
#'                                       archiMultipartite, connectParam, model,
#'                                       directed, dimLabels = c('A','B'), seed = 2)
#' listSBM <- mySampleMSBM$listSBM
#' estimOptions <- list(initBM = FALSE, nbCores  = 2)
#' myMSBM <- estimateMultipartiteSBM(listSBM, estimOptions)
#' plot(myMSBM, type = "data")
#' plot(myMSBM, type = "expected")
#' plot(myMSBM, type = "meso")
#' }
estimateMultipartiteSBM <- function(listSBM,
                                    estimOptions = list())
{

  myMSBM <- MultipartiteSBM_fit$new(listSBM)


  currentOptions <- list(
    verbosity     = 1,
    nbBlocksRange = rep(list(c(1,10)), length(myMSBM$dimLabels)),
    nbCores       = 2,
    maxiterVE     = 100,
    maxiterVEM    = 100,
    initBM = TRUE
  )

  names(currentOptions$nbBlocksRange) <- myMSBM$dimLabels
  ## Current options are default expect for those passed by the user
  currentOptions[names(estimOptions)] <- estimOptions

  myMSBM$optimize(currentOptions)

  myMSBM
}




# -------------------------------------------------------------------------
#' Estimation for Multiplex SBM
#'
#' @param listSBM list of networks that were defined by the \code{defineSBM} function
#' @param dependent logical parameter indicating whether the networks in the multiplex structure are dependent beyond the latent variables,
#' @param estimOptions options for the inference procedure
#' @details The list of parameters \code{estimOptions} essentially tunes the optimization process and the variational EM algorithm, with the following parameters
#' * "nbCores" integer for number of cores used.  Default is 2
#' * "verbosity" integer for verbosity (0, 1). Default is 1
#' * "nbBlocksRange" List of length the number of functional groups, each element supplying the minimal and maximal number of blocks to be explored. The names of the list must be the names of the functional groups.  Default value is from 1 to 10)
#' * "initBM" Boolean. True if using simple and bipartite SBM as initialisations. Default value  = TRUE
#' * "maxiterVEM" Number of max. number of iterations in  the VEM. Default value  = 100
#' * "maxiterVE" Number of max. number of iterations in  the VE. Default value  = 100
#' * "plot" boolean, should the ICL by dynamically plotted or not. Default is TRUE. For dependent networks
#' * "exploreFactor" control the exploration of the number of groups. For dependent networks
#' * "exploreMin" explore at least until exploreMin even if the exploration factor rule is achieved. Default 4. See the package blockmodels for details. For dependent networks
#' * "exploreMax" Stop exploration at exploreMax  even if the exploration factor rule is not achieved. Default Inf. See the package blockmodels for details. For dependent networks
#' * "nbBlocksRange" minimal and maximal number or blocks explored. For dependent networks
#' * "fast" logical: should approximation be used for Bernoulli model with covariates. Default to \code{TRUE}. For dependent networks
#'
#' @return a MultiplexSBM_fit object with the estimated parameters and the blocks
#' @export
#'
#' @examples
#' \dontrun{
#' ### =======================================
#' ### MULTIPLEX SBM without dependence between layers
#' ##
#' Nnodes <- 40
#' blockProp <- c(.4,.6)
#' nbLayers <- 2
#' connectParam <- list(list(mean=matrix(rbeta(4,.5,.5),2,2)),list(mean=matrix(rexp(4,.5),2,2)))
#' model <- c("bernoulli","poisson")
#' type <- "directed"
#' mySampleMultiplexSBM <-
#'    sampleMultiplexSBM(
#'    nbNodes = Nnodes,
#'     blockProp = blockProp,
#'    nbLayers = nbLayers,
#'    connectParam = connectParam,
#'    model=model,
#'    type=type)
#' listSBM <- mySampleMultiplexSBM$listSBM
#' estimOptions <- list(initBM = FALSE, nbCores  = 2)
#' myMultiplexSBM <- estimateMultiplexSBM(listSBM,estimOptions,dependent=FALSE)

#' ### =======================================
#' ### MULTIPLEX SBM Gaussian with dependence
#' ##
#' Q <- 3
#' nbLayers <- 2
#' connectParam <- list()
#' connectParam$mu <- vector("list",nbLayers)
#' connectParam$mu[[1]] <-  matrix(.1,Q,Q) + diag(1:Q)
#' connectParam$mu[[2]] <- matrix(-2,Q,Q) + diag(rev(Q:1))
#' connectParam$Sigma <- matrix(c(2,1,1,4),nbLayers,nbLayers)
#' model <- rep("gaussian",2)
#' type <- "directed"
#' Nnodes <- 80
#' blockProp <- c(.3,.3,.4)
#' mySampleMultiplexSBM <-
#'   sampleMultiplexSBM(
#'      nbNodes = Nnodes,
#'      blockProp = blockProp,
#'      nbLayers = nbLayers,
#'      connectParam = connectParam,
#'      model=model,
#'      type="undirected",
#'      dependent=TRUE)
#' listSBM <- mySampleMultiplexSBM$listSBM
#' myMultiplexSBM <- estimateMultiplexSBM(listSBM,estimOptions,dependent=TRUE)
#' ## MultiplexSBM Bernoulli with dependence
#' Q <- 2
#' P00<-matrix(runif(Q*Q),Q,Q)
#' P10<-matrix(runif(Q*Q),Q,Q)
#' P01<-matrix(runif(Q*Q),Q,Q)
#' P11<-matrix(runif(Q*Q),Q,Q)
#' SumP<-P00+P10+P01+P11
#' P00<-P00/SumP
#' P01<-P01/SumP
#' P10<-P10/SumP
#' P11<-P11/SumP
#' connectParam <- list()
#' connectParam$prob00 <- P00
#' connectParam$prob01 <- P01
#' connectParam$prob10 <- P10
#' connectParam$prob11 <- P11
#' model <- rep("bernoulli",2)
#' type <- "directed"
#' nbLayers <- 2
#' Nnodes <- 40
#' blockProp <- c(.6,.4)
#' mySampleMultiplexSBM <-
#'    sampleMultiplexSBM(
#'      nbNodes = Nnodes,
#'      blockProp = blockProp,
#'      nbLayers = nbLayers,
#'      connectParam = connectParam,
#'      model=model,
#'      type=type,
#'      dependent=TRUE)
#' listSBM <- mySampleMultiplexSBM$listSBM
#' myMultiplexSBM <- estimateMultiplexSBM(listSBM,estimOptions,dependent=TRUE)
#' }
estimateMultiplexSBM <- function(listSBM,
                                    dependent = FALSE,
                                    estimOptions = list())
{


  myMSBM <- MultiplexSBM_fit$new(listSBM,dependentNet= dependent)


  if (dependent)
  {
    if (!is.null(estimOptions$verbosity)){
      if ((estimOptions$verbosity == 0)&(is.null(estimOptions$plot))){estimOptions$plot   = FALSE}
    }
    currentOptions <- list(
      verbosity     = 3,
      plot          = TRUE,
      explorFactor  = 1.5,
      nbBlocksRange = c(4,Inf),
      nbCores       = 2,
      fast          = TRUE
    )
  } else
  {
    currentOptions <- list(
      verbosity     = 1,
      nbBlocksRange = rep(list(c(1, 10)), length(myMSBM$dimLabels)),
      nbCores       = 2,
      maxiterVE     = 100,
      maxiterVEM    = 100,
      initBM = TRUE
    )
  }
  names(currentOptions$nbBlocksRange) <- myMSBM$dimLabels
  ## Current options are default expect for those passed by the user
  currentOptions[names(estimOptions)] <- estimOptions

  myMSBM$optimize(currentOptions)

  myMSBM

}
