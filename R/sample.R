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
#' @return  an object with class \code{\link{SimpleSBM}}
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
  mySampler$rNetwork(store = TRUE)
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
#' @return an object with class \code{\link{BipartiteSBM}}
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
#' mySampler$rEdges()   # sample new edges
#' mySampler$rNetwork()   # sample a new networrk (blocks and edges)
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
  mySampler$rNetwork(store = TRUE)
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
#' @return  a list of two elements : \code{simulatedMemberships} are the clustering of each node in each Functional Group,
#'   \code{multipartiteNetwork} is the list of the simulated networks (each one being  a simple or bipartite network)
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


 # browser()
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

#' Sampling of Multiplex SBMs
#'
#' This function samples a Multiplex Stochastic Block Models, with various model
#' for the distribution of the edges:  Bernoulli, Poisson, or Gaussian models
#'
#' @param nbNodes number of nodes in each functional group involved in the Multiplex network
#' @param blockProp a vector for block proportion if the networks are simple, a list of parameters for block proportions for both functional groups if the networks are bipartite
#' @param nbLayers a matrix with two columns and nbNetworks lines, each line specifying the index of the functional groups in interaction.
#' @param connectParam list of parameters for connectivity (of length nbNetworks). Each element is a list of one or two elements: a matrix of means 'mean' and an optional matrix of variances 'var', the sizes of which must match \code{blockProp} length
#' @param model a vector of characters describing the model for  each network of the Multiplex relation between nodes (\code{'bernoulli'}, \code{'poisson'}, \code{'gaussian'}, ...). Default is \code{'bernoulli'}.
#' @param type a string of character indicating whether the networks are directed, undirected or bipartite
#' @param dependent connection parameters in each network
#' @param dimLabels an optional list of labels for functional group involved in the network
#' @param seed numeric to set the seed.
#' @return  a list of two elements : \code{simulatedMemberships} are the clustering of each node in each Functional Group,  \code{MultiplexNetwork} is the list of the simulated networks (each one being  a simple or bipartite network)
#'
#' @examples
#' nbLayers <- 2
#'
#' ## MultiplexSBM without dependence between layers
#' Nnodes <- 40
#' blockProp <- c(.4,.6)
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
#'
#' ## MultiplexSBM Gaussian with dependence
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
#' listSBM_BB <- mySampleMultiplexSBM$listSBM
#'
#' @importFrom stats rmultinom runif rnorm
#' @export
sampleMultiplexSBM <- function(nbNodes,
                               blockProp,
                               nbLayers,
                               connectParam,
                               model,
                               type=c("directed","undirected","bipartite"),
                               dependent=FALSE,
                               dimLabels = NULL,
                               seed = NULL) {


  # dimLabels for compatibility with define SBM
  if (is.null(dimLabels))
  {
    if (type=="bipartite") dimLabels = c(row = "row", col = "col")
    else dimLabels = "node"
  }



  # block prop list ou simple vecteur
  if (!((length(nbNodes)==2 & is.list(blockProp)) | (length(nbNodes)==1 & !is.list(blockProp)) | (length(nbNodes)==1 & is.list(blockProp) & length(blockProp)==1 )))
    stop("length of vector nbNodes and length of list blockProp should match")

  # same sanity check as in the R6 class MultiplexSBM_fit
  # check whether the Multiplex at hand is actually a multiplex
  # if (any(c(unique(archiMultiplex[,1]),unique(archiMultiplex[,2])) > 1))
  #   stop("Architecture of networks provided does not correspond to a Multiplex architecture")

  # CHECKING dependence structure
  if (dependent) {
    # on empeche cette option
    # if (! ( all(directed == TRUE) | all(directed == FALSE)) )
    #   stop("in the dependent case, all networks should be either directed or not directed")

    dBern  <- isTRUE(all.equal(model, rep("bernoulli", length(model))))
    dGauss <- isTRUE(all.equal(model, rep("gaussian" , length(model))))
    if (!(dGauss | (dBern&length(model) == 2)))
      stop("dependency in multiplex network is only handled for Gaussian distribution or a bivariate Bernoulli distribution")

    if (dBern)
    {
      P00 <- connectParam$prob00
      P01 <- connectParam$prob01
      P10 <- connectParam$prob10
      P11 <- connectParam$prob11

      if (type == "bipartite")
      {
        Z1 <- t(rmultinom(nbNodes[1], size = 1, prob = blockProp[[1]]))
        Z2 <- t(rmultinom(nbNodes[2], size = 1, prob = blockProp[[2]]))
        MU <-matrix(runif(prod(nbNodes)),nbNodes[1],nbNodes[2])
        M1 <-1*(MU>Z1%*%(P00+P01)%*%t(Z2))
        M2 <-1*(((MU>Z1%*%(P00)%*%t(Z2)) & (MU<Z1%*%(P00+P01)%*%t(Z2))) | (MU>Z1%*%(1-P11)%*%t(Z2)))
        memberships <- list(as_clustering(Z1),as_clustering(Z2))
      }
      else {
        if (!is.list(blockProp)) blockProp = list(blockProp)
        Z <- t(rmultinom(nbNodes[1], size = 1, prob = blockProp[[1]]))
        MU <-matrix(runif((nbNodes)**2),nbNodes,nbNodes)
        M1 <-1*(MU>Z%*%(P00+P01)%*%t(Z))
        M2 <-1*(((MU>Z%*%(P00)%*%t(Z)) & (MU<Z%*%(P00+P01)%*%t(Z))) | (MU>Z%*%(1-P11)%*%t(Z)))
        memberships <- list(as_clustering(Z))
      if (type== "undirected")
      {
        if (any(!sapply(connectParam,isSymmetric))) stop("Non symmetric parameters")
        MU[lower.tri(MU)]<-t(MU)[lower.tri(MU)]
        M1 <-1*(MU>Z%*%(P00+P01)%*%t(Z))
        M2 <-1*(((MU>Z%*%(P00)%*%t(Z)) & (MU<Z%*%(P00+P01)%*%t(Z))) | (MU>Z%*%(1-P11)%*%t(Z)))
      }
      }

      listNetworks <- list()
      names(memberships) <- dimLabels
      type_l <- ifelse(type=="bipartite","bipartite","simple")
      listNetworks[[1]] <- defineSBM(netMat  = M1, model = model[1], type = type_l,dimLabels =  dimLabels)
      listNetworks[[2]] <- defineSBM(netMat  = M2, model = model[2], type = type_l,dimLabels =  dimLabels)

    }
    if (dGauss)
    {
      listNetworks <- vector("list",nbLayers)
      Mus <- connectParam$mu
      Sig <- connectParam$Sigma

      if (type == "bipartite")
      {
        Z1 <- t(rmultinom(nbNodes[1], size = 1, prob = blockProp[[1]]))
        Z2 <- t(rmultinom(nbNodes[2], size = 1, prob = blockProp[[2]]))
        Noise <- t(chol(Sig)) %*% matrix(rnorm(prod(nbNodes)*nbLayers),nbLayers,prod(nbNodes))

        for (l in 1:nbLayers)
        {
            nettemp <- Z1%*%Mus[[l]]%*%t(Z2) + matrix(Noise[l,],nbNodes[1],nbNodes[2])
            listNetworks[[l]] <- defineSBM(netMat  = nettemp, model = model[l], type = "bipartite",dimLabels =  dimLabels)
        }
        memberships <- list(as_clustering(Z1),as_clustering(Z2))
      }
      else
      {
        if (!is.list(blockProp)) blockProp = list(blockProp)
        Z <- t(rmultinom(nbNodes[1], size = 1, prob = blockProp[[1]]))
        Noise <- t(chol(Sig)) %*% matrix(rnorm(nbNodes*nbNodes*nbLayers),nbLayers,nbNodes*nbNodes)

        for (l in 1:nbLayers)
        {
          nettemp <- Z%*%Mus[[l]]%*%t(Z) + matrix(Noise[l,],nbNodes,nbNodes)
          if (type== "undirected")
          {
            nettemp[lower.tri(nettemp)] <- t(nettemp)[lower.tri(nettemp)]
          }
          listNetworks[[l]] <- defineSBM(netMat  = nettemp, model = model[l], type = "simple",dimLabels =  dimLabels)
        }
        memberships <- list(as_clustering(Z))
      }
      names(memberships) <- dimLabels

    }

    return(  list(listSBM =  listNetworks, memberships  = memberships))
    }


   if (!dependent)  {
     if (length(unique(c(length(model),length(connectParam),nbLayers)))>1)
       stop("length of vector model and length of list connectParam and number of layers should match")

     archiMultipartite <- matrix(1,nrow=nbLayers,ncol=2)
     if (type == "bipartite") archiMultipartite[,2] <- 2
     directed <- rep(type=="directed",nbLayers)
     if (!is.list(blockProp)) blockProp = list(blockProp)
     print("use of sampleMultipartite")
     return(sampleMultipartiteSBM(nbNodes,blockProp,archiMultipartite,connectParam,model,directed,dimLabels,seed))
   }

}
