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
#' @param dimLabels an optional list of labels for each dimension (in row, in column)
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
#' plot(mySampler,type='meso')
#' hist(mySampler$networkData)
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
#' plot(mySampler,type='meso')
#' hist(mySampler$networkData)
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
#' plot(mySampler,type='meso')
#' hist(mySampler$networkData)
#' @export
sampleSimpleSBM <- function(nbNodes,
                            blockProp,
                            connectParam,
                            model = 'bernoulli',
                            directed = FALSE,
                            dimLabels = c(node = "nodeName"),
                            covariates = list(),
                            covariatesParam = numeric(0)) {

  mySampler <- SimpleSBM$new(model, nbNodes, directed, blockProp, connectParam, dimLabels, covariatesParam, covariates)
  mySampler$rMemberships(store = TRUE)
  mySampler$rAdjacency(store = TRUE)
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
#' @param model character describing the model for the relation between nodes (\code{'bernoulli'}, \code{'poisson'}, \code{'gaussian'}, \code{'ZIgaussian'}). Default is \code{'bernoulli'}.
#' @param dimLabels an optional list of labels for each dimension (in row, in column)
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
#' dimLabels = c(row='Reader',col='Book')
#' mySampler <- sampleBipartiteSBM(nbNodes, blockProp, connectParam, model = 'bernoulli',dimLabels)
#' plot(mySampler)
#' plot(mySampler,type='meso',plotOptions = list(vertex.label.name=list(row='Reader',col='Book')))
#' plot(mySampler,type='meso',plotOptions = list(vertex.label.name=c('A','B'),vertex.size = 1.4))
#' mySampler$rMemberships() # sample new memberships
#' mySampler$rIncidence()   # sample new incidence matrix

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
#' dimLabels = c(row = 'Ind', col = 'Service')
#' mySampler <- sampleBipartiteSBM(nbNodes, blockProp, connectParam, model = 'poisson', dimLabels)
#' plot(mySampler,type='expected')
#' plotOptions = list(vertex.label.name=c('U','V'),vertex.size = c(1.4,1.3))
#' plot(mySampler, type='meso', plotOptions = plotOptions)
#' hist(mySampler$networkData)
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
#' hist(mySampler$networkData)
#'
#' @export
sampleBipartiteSBM <- function(nbNodes,
                            blockProp,
                            connectParam,
                            model = 'bernoulli',
                            dimLabels    = c(row = "rowName", col = "colName"),
                            covariates = list(),
                            covariatesParam = numeric(0)) {

  mySampler <- BipartiteSBM$new(model, nbNodes, blockProp, connectParam, dimLabels, covariatesParam, covariates)
  mySampler$rMemberships(store = TRUE)
  mySampler$rIncidence(store = TRUE)
  mySampler
}
#' Sampling of Multipartite SBMs
#'
#' This function samples a Multipartite Stochastic Block Models, with various model
#' for the distribution of the edges:  Bernoulli, Poisson, or Gaussian models
#'
#' @param nbNodes number of nodes in each functional group involved in the multipartite network
#' @param blockProp a list of parameters for block proportions  in each functional group
#' @param archiMultipartite a matrix with two columns and nbNetworks lines, each line specifying the index of the functional groups in interaction.
#' @param connectParam list of parameters for connectivity (of length nbNetworks). Each element is a list of one or two elements: a matrix of means 'mean' and an optional matrix of variances 'var', the sizes of which must match \code{blockProp} length
#' @param model a vector of characters describing the model for  each network of the Multipartite relation between nodes (\code{'bernoulli'}, \code{'poisson'}, \code{'gaussian'}, ...). Default is \code{'bernoulli'}.
#' @param directed a vector of logical, directed network or not for each network. Default is \code{FALSE}.
#' @param dimLabels an optional list of labels for functional group involved in the network
#' @param seed numeric to set the seed.
#' @return  a list of two elements : \code{simulatedMemberships} are the clustering of each node in each Functional Group,  \code{multipartiteNetwork} is the list of the simulated networks (each one being  a simple or bipartite network)
#'
#' @examples
#' ### =======================================
#' ### MULTIPARTITE SBM  : 4 networks between 3 Functional Groups
#' ## Graph parameters
#' # About the Functional Groups (FG)
#' nbNodes <-  c(100,50,40)
#' blockProp <- vector("list", 3)  # parameters of clustering in each functional group
#' blockProp[[1]] <- c(0.4,0.3,0.3) # in Functional Group 1
#' blockProp[[2]] <- c(0.6,0.4) # in Functional Group 2
#' blockProp[[3]]  <- c(0.6,0.4) # in Functional Group 3
#' # About the interactions between the FG
#' archiMultipartite  <-  rbind(c(1,2),c(2,3),c(2,2),c(1,3)) #
#' model <- c('bernoulli','poisson','gaussian','gaussian') # type of distribution in each network
#' # for each network : directed or not (not required for an interaction between two different FG)
#' directed <- c( NA, NA  ,  FALSE , NA)
#' connectParam <- list()
#' connectParam[[1]] <- list(mean = matrix(c(0.3, 0.3, 0.5, 0.2, 0.6, 0.6),3,2))
#' connectParam[[2]] <- list(mean = matrix(c(1000 , 500,  400 , 950),2,2))
#' connectParam[[3]] <- list(mean = matrix(c(10, 0, -10, 20), 2,2), var = matrix(1,2,2))
#' connectParam[[4]] <- list(mean = matrix(c(3, 23 ,11 ,16 , 2 ,25), 3,2))
#' connectParam[[4]]$var <- matrix(c(10,20,1,5,0.1,10), 3,2)
#' dimLabels <- c('A','B','C')
#' ## Graph Sampling
#' mySampleMBM <- sampleMultipartiteSBM(nbNodes, blockProp,
#'                                      archiMultipartite,
#'                                      connectParam, model, directed,
#'                                      dimLabels,seed = 3)
#' listSBM <- mySampleMBM$listSBM
#' memberships <- mySampleMBM$memberships
#' plotMyMultipartiteMatrix(listSBM)
#' plotMyMultipartiteMatrix(listSBM,plotOptions = list(normalized = TRUE))
#' plotMyMultipartiteMatrix(listSBM,memberships = memberships,plotOptions = list(normalized = TRUE))
#' @export
sampleMultipartiteSBM <- function(nbNodes,
                            blockProp,
                            archiMultipartite,
                            connectParam,
                            model,
                            directed,
                            dimLabels = NULL,
                            seed = NULL) {

  nbNetworks <- length(connectParam)

  # transfo intro GREMLINS param
  list_theta <- list()
  for (l in 1 : nbNetworks){
    if  ( model[l] %in% c('bernoulli','poisson')) { list_theta[[l]] = connectParam[[l]]$mean} else{list_theta[[l]] =  connectParam[[l]]}
  }
  list_pi  <- blockProp
  v_NQ <- nbNodes
  E <- archiMultipartite
  typeInter <- rep(0,nbNetworks)
  for (l in 1 : nbNetworks){
    if (E[l,2] != E[l,1]){typeInter[l] = 'inc'} else {
      if (directed[l]){typeInter[l] = 'diradj'} else{typeInter[l]  = 'adj'}
    }
  }
  v_distrib <- model
  namesFG <-  dimLabels
  dataSimGREMLIN <- rMBM(v_NQ ,E,typeInter,  v_distrib, list_pi , list_theta, namesFG, keepClassif = TRUE, seed= seed)
  listNetworks <- list()
  memberships <- dataSimGREMLIN$classif
  names(memberships) <- namesFG

  for (l in 1:nbNetworks){
    dimLabels_l = c(row = dataSimGREMLIN$list_Net[[l]]$rowFG, col = dataSimGREMLIN$list_Net[[l]]$colFG)
    type_l <- ifelse(typeInter[l] == 'inc','bipartite','simple')
    if (type_l == "simple") dimLabels_l = c(node = unique(dimLabels_l))
    listNetworks[[l]] <- defineSBM(netMat  = dataSimGREMLIN$list_Net[[l]]$mat, model = model[l], type = type_l, directed = directed[l],dimLabels =  dimLabels_l)
  }

  list(listSBM =  listNetworks, memberships  = memberships)

}

