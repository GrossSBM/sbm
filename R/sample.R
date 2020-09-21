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
#' plot(mySampler,type='meso')
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
#' plot(mySampler,type='meso')
#' hist(mySampler$netMatrix)
#' @export
sampleSimpleSBM <- function(nbNodes,
                            blockProp,
                            connectParam,
                            model = 'bernoulli',
                            directed = FALSE,
                            dimLabels    = list(row = "rowLabel", col = "colLabel"),
                            covariates = list(),
                            covariatesParam = numeric(0)) {

  mySampler <- SimpleSBM_sampler$new(model, nbNodes, directed, blockProp, connectParam, dimLabels, covariatesParam, covariates)
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
#' mySampler <- sampleBipartiteSBM(nbNodes, blockProp, connectParam, model = 'bernoulli',dimLabels=c('Readers','Book'))
#' plot(mySampler)
#' mySampler$rMemberships() # sample new memberships
#' mySampler$rIncidence()   # sample new incidence matrix
#' plot(mySampler,type='meso',plotOptions=list(vertex.size=1.4))
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
                            dimLabels    = list(row = "rowLabel", col = "colLabel"),
                            covariates = list(),
                            covariatesParam = numeric(0)) {

  mySampler <- BipartiteSBM_sampler$new(model, nbNodes, blockProp, connectParam, dimLabels, covariatesParam, covariates)
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
#'
#' @return  a list of two elements : \code{simulatedMemberships} are the clustering of each node in each Functional Group,  \code{multipartiteNetwork} is the list of the simulated networks (each one being  a simple or bipartite network)
#'
#' @examples
#' ### =======================================
#' ### MULTIPARTITE SBM  : 4 networks between 3 Functional Groups
#' ## Graph parameters
#' # About the Functional Groups (FG)
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
#' connectParam[[3]] <- list(mean  =  matrix(rbeta(nbBlocks[E[3,1]] * nbBlocks[E[3,2]],0.9,0.9 ), nrow = nbBlocks[E[3,1]], ncol = nbBlocks[E[3,2]]))
#' connectParam[[3]]$mean <-  0.5*(connectParam[[3]]$mean + t(connectParam[[3]]$mean)) # symetrisation for network 3
#' connectParam[[4]] <- list(mean = matrix(rnorm(nbBlocks[E[4,1]] * nbBlocks[E[4,2]],7.5,10 ), nrow = nbBlocks[E[4,1]], ncol = nbBlocks[E[4,2]]))
#' connectParam[[4]]$var <- matrix(rgamma(nbBlocks[E[4,1]] * nbBlocks[E[4,2]],7.5,0.1 ), nrow = nbBlocks[E[4,1]], ncol = nbBlocks[E[4,2]])
#' dimLabels <- as.list(c('A','B','C'))
#' seed <- 10
#' ## Graph Sampling
#' mySampleMBM <- sampleMultipartiteSBM(nbNodes, blockProp, archiMultipartite, connectParam, model, directed, dimLabels,seed)
#' listSBM <- mySampleMBM$listSBM
#' memberships <- mySampleMBM$memberships
#' plotMyMultipartiteMatrix(listSBM,memberships)
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

  # transfo intro GREMLIN param
  list_theta <- list()
  for (l in 1 : nbNetworks){
    if  ( model[l] != 'gaussian') { list_theta[[l]] = connectParam[[l]]$mean} else{list_theta[[l]] =  connectParam[[l]]}
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
  namesFG <-  unlist(dimLabels)
  dataSimGREMLIN <- rMBM(v_NQ ,E,typeInter,  v_distrib, list_pi , list_theta, namesFG, keepClassif = TRUE, seed= seed)
  listNetworks <- list()
  memberships <- dataSimGREMLIN$classif
  names(memberships) <- namesFG
  for (l in 1:nbNetworks){
    dimLabels_l = list(row = dataSimGREMLIN$list_Net[[l]]$rowFG, col = dataSimGREMLIN$list_Net[[l]]$colFG)
    type_l <- ifelse(typeInter[l] == 'inc','bipartite','simple')
    listNetworks[[l]] <- defineSBM(netMat  = dataSimGREMLIN$list_Net[[l]]$mat, model = model[l], type = type_l, directed = directed[l],dimLabels =  dimLabels_l)
  }

  list(listSBM =  listNetworks, memberships  = memberships)



}

