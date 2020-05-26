############################################################
############ Basic functions to use in ScoreSBM ###############


#-------------------------------------------------------------------------------
# Computes the number of (symmetric) pairs
#-------------------------------------------------------------------------------
n2nbPairs <- function(n, symmetric, diag=FALSE){
  if(symmetric){
    N <- ifelse(diag, n*(n+1)/2, n*(n-1)/2)
  }else{
    N <- ifelse(diag, N <- n^2, n*(n-1))
  }
  return(N)
}

#-------------------------------------------------------------------------------
# Computes the individual from the number of (symmetric) pairs
#-------------------------------------------------------------------------------
nbPairs2n <- function(N, symmetric, diag=FALSE){
  if(symmetric){
    n <- ifelse(diag, (-1 + sqrt(1 + 8*N))/2 , (1 + sqrt(1 + 8 * N))/2)
  }else{
    n <- ifelse(diag, sqrt(N) , (1 + sqrt(1 + 4*N))/2)
  }
  return(n)
}

#-------------------------------------------------------------------------------
# Builds a (symmetric) matrix from a vector
#-------------------------------------------------------------------------------
vect2Mat <- function(V, symmetric, diag=FALSE)
{

  N <- length(V); n <- nbPairs2n(N, symmetric=symmetric, diag=diag);
  M <- matrix(0,n,n)

  if (symmetric == TRUE){
    # n <- ifelse(diag, (-1 + sqrt(1 + 8*N))/2 , (1 + sqrt(1 + 8 * N))/2)
    # M <- matrix(0,n,n)
    M[lower.tri(M, diag = diag)] <- V
    B <- t(M)
    diag(B) <- 0
    M <- B + M
  }else{
    # n <- ifelse(diag, sqrt(N) , (1 + sqrt(1 + 4*N))/2)
    # M <- matrix(0,n,n)
    where <- (lower.tri(M, diag = diag) +  upper.tri(M, diag = diag)) > 0
    M[where] <- V
  }
  return(M)
}

#-------------------------------------------------------------------------------
# Builds a vector from a (symmetric) matrix
#-------------------------------------------------------------------------------
mat2Vect <- function(M, symmetric, diag=FALSE)
{
  if (symmetric) {
    where <- lower.tri(M, diag = diag)
  } else {
    where <- (lower.tri(M, diag=diag) +  upper.tri(M, diag = diag)) > 0
  }
  V <- M[where]
  return(V)
}

#-------------------------------------------------------------------------------
#- Gest the indices corresponding to a list of (symmetric) pairs
#-------------------------------------------------------------------------------
indices <- function(n, symmetric, diag=FALSE)
{
  N <- n2nbPairs(n, symmetric=symmetric, diag=diag)
  S <- vect2Mat(1:N, symmetric=symmetric, diag=diag)
  if (symmetric) {S[upper.tri(S)] = 0}
  res <- which(S!=0, arr.ind=TRUE)
  return(res)
}



#-------------------------------------------------------------------------------
#---- Transorm the list of the matrices into a unique matrix
#-------------------------------------------------------------------------------
scoreList2scoreMat <- function(listScores,symmetric){

  S <- listScores
  d <- length(S);
  #n <- nrow(S[[1]])
  mat_S <- sapply(1:d, function(q) {mat2Vect(S[[q]], symmetric = symmetric, diag = F)})
  # of dimension N  * d where N = n(n-1)/2 if symmetric, n(n-1) otherwise
  return(mat_S)
}


#-------------------------------------------------------------------------------
#---------- Distance on tau--------------------------------------
#-------------------------------------------------------------------------------

distTau  <- function(tau,tauOld)
{
  return(sqrt(sum(as.vector(tau - tauOld)^2)))
}


#------------------------------------------------------------------------------
#' Initialization of the inference procedure
#-------------------------------------------------------------------------------
#' This function initialises the inference method by mixing a Gaussian mixture on the scores of each pair of nodes and a SBM on the resulting estimated network G
#' @param scoreList a list of the Scores (matrices of size nbNodes x nbNodes)
#' @param directed  a logical : TRUE if the underlying network is directed,  FALSE otherwise (default value FALSE).
#' @param estimOptions a list of parameters controlling the initialisation step of the inference method. See details.
#'
initInferenceScoreSBM <- function(scoreList, directed = FALSE, estimOptions = list()){

  currentOptions <- list(
    explorFactor  = 1.5,
    nbBlocksRange = c(4,Inf),
    nbCores       = parallel::detectCores()
  )

  ## Current options are default expect for those passed by the user
  currentOptions[names(estimOptions)] <- estimOptions



  #------------------------------------------------
  nbScores <- length(scoreList)
  scoreMat <- sapply(1:nbScores , function(q) {mat2Vect(scoreList[[q]], symmetric = !directed, diag = F)})

  N <- nrow(scoreMat)
  #--------  Mclust on the d scores dyad per dyad
  param_gm <- mclust::Mclust(scoreMat, G = 2, verbose = FALSE)
  psi <- param_gm$z
  G <- param_gm$classification - 1
  mu  <- param_gm$parameters$mean
  test_G <- rowMeans(t(mu)) #identify G = 0  and G =1
  if (test_G[1] > test_G[2]) {
    psi <- psi[,c(2,1)]
    G = 1 - G
  }


  #------------------ init of SBM parameters
  membership_type <- ifelse(directed, "SBM", "SBM_sym")
  param_sbm <- blockmodels::BM_bernoulli(membership_type, adj = vect2Mat(G, symmetric = !directed),
                                         plotting = '',
                                         verbosity = 0,
                                         explore_min = currentOptions$nbBlocksRange[1],
                                         explore_max = currentOptions$nbBlocksRange[2],
                                         exploration_factor = currentOptions$explorFact,
                                         ncores = currentOptions$nbCores)

  param_sbm$estimate()
  Kmax <- length(param_sbm$memberships)
  initAll <- initInferenceScoreSBM(scoreList, directed,currentOptions)

  #-------------------------------------------------------------------------
  #--------------------ESTIMATION ALL MODELS -------------------------------
  #------------------------------------------------------------------------
  nbModels <- length(initAll$tau)
  nbBlocksAll = vapply(1:nbModels,function(m){ncol(initAll$tau[[m]])},1)
  indModels <- (1:nbModels)[nbBlocksAll <= currentOptions$nbBlocksRange[2]]
  indModels <- (1:nbModels)[nbBlocksAll >= currentOptions$nbBlocksRange[1]]
  nbBlocksAll <- nbBlocksAll[indModels]


  tau_init <-  lapply(1:Kmax, function(K){param_sbm$memberships[[K]]$Z})
  eta  <-  lapply(1:Kmax,  function(K){array(rep(psi[, 2], K * K), c(N, K, K))})
  #-------------------------------------------------------------------------------

  res <- list(psi = psi, tau = tau_init, eta = eta, ICL = param_sbm$ICL, G = G)
  return(res)

}

#-------------------------------------------------------------------------------
# M step of the VEM algorithm
#-------------------------------------------------------------------------------
#' @importFrom mvtnorm dmvnorm
mStepScoreSBM <- function(scoreMat, qDist, directed){

  # scoreMat <- scoreMat; qDist <- qDist; directed <- FALSE

  # Dimensions
  d <- ncol(scoreMat); nbBlocks <- ncol(qDist$tau) #N <- nrow(scoreMat);

  # Proportions
  blockProp <- colMeans(qDist$tau)

  # Connection probabilities
  connectParam <- matrix(0, nbBlocks, nbBlocks)
  sapply(1:nbBlocks, function(k){sapply(1:nbBlocks, function(l){
    tauVec <- mat2Vect(qDist$tau[, k] %o% qDist$tau[, l], symmetric = !directed, diag = FALSE)
    connectParam[k, l] <<- tauVec %*% qDist$eta[, k, l] / sum(tauVec)
  })})

  # Emission distributions: mu and Sigma
  mu <- matrix(0, 2, d);
  Sigma <- array(dim = c(2, d, d))
  sapply(1:2, function(g){
    mu[g, ] <<- t(qDist$psi[, g]) %*% scoreMat / sum(qDist$psi[, g])
    Sigma[g, , ] <<- t(scoreMat) %*% diag(qDist$psi[, g]) %*% scoreMat / sum(qDist$psi[, g])
    Sigma[g, , ] <<- Sigma[g, , ] - mu[g, ] %o% mu[g, ]
    Sigma[g, , ] <<- .5 * (Sigma[g, , ] + t(Sigma[g, , ]))
  })
  emissionParam <- list(noEdgeParam = list(mean = mu[1, ], var = Sigma[1, , ]),
                        EdgeParam = list(mean = mu[2, ], var = Sigma[2, , ]))

  res <- list(blockProp = blockProp, connectParam = connectParam, emissionParam = emissionParam)
  return(res)
}

#-------------------------------------------------------------------------------
# VE step of the VEM algorithm
#-------------------------------------------------------------------------------
veStepScoreSBM <- function(scoreMat, theta,tauOld, directed, estimOptions = list()){

  # scoreMat <- scoreMat; theta <- thetaHat; directed <- FALSE
  # epsilon_tau <- epsilon_eta <- 1e-4; tauOld <- qDist$tau
  currentOptions <- list(
    maxIterVE = 100 ,
    tauTol = 2 * .Machine$double.eps,
    valStopCrit = 1e-6,
    etaTol = 2 * .Machine$double.eps
  )
  currentOptions[names(estimOptions)] <- estimOptions

  noConvergence = 0
  # Dimensions
  nbBlocks <- length(theta$blockProp);
  N <- nrow(scoreMat); n <- nbPairs2n(N, symmetric = !directed)
  indexList <- indices(n, symmetric = !directed)


  # log(phi)
  logPhi <- matrix(0, N, 2)
  logPhi[, 1] <- mvtnorm::dmvnorm(scoreMat,
                                  mean = theta$emissionParam$noEdgeParam$mean,
                                  sigma = theta$emissionParam$noEdgeParam$var, log = TRUE)
  logPhi[, 2] <- mvtnorm::dmvnorm(scoreMat,
                                  mean = theta$emissionParam$EdgeParam$mean,
                                  sigma = theta$emissionParam$EdgeParam$var, log = TRUE)

  # eta
  eta <- array(dim = c(N, nbBlocks, nbBlocks))
  invisible(sapply(1:nbBlocks, function(k){sapply(1:nbBlocks, function(l){ # k <- 1; l <- 2
    etaTmp <- logPhi + (rep(1, N) %o% c(log(1 - theta$connectParam[k, l]), log(theta$connectParam[k, l])))
    etaTmp <- etaTmp - apply(etaTmp, 1, max)
    etaTmp <- exp(etaTmp); etaTmp <- etaTmp / rowSums(etaTmp)
    etaTmp <- etaTmp + currentOptions$etaTol; etaTmp <- etaTmp / rowSums(etaTmp)
    eta[, k, l] <<- etaTmp[, 2]
  })}))

  # log(A)
  logA <- array(dim = c(N, nbBlocks, nbBlocks))
  sapply(1:nbBlocks, function(k){sapply(1:nbBlocks, function(l){ # k <- 1; l <- 2
    logA[, k, l] <<- (1 - eta[, k, l])*(log(1 - theta$connectParam[k, l]) + logPhi[, 1] - log(1-eta[,k,l])) +
      eta[, k, l]*(log(theta$connectParam[k, l]) + logPhi[, 2] -  log(eta[,k,l]))
  })})


  #-------------- Fixed point
  tau <- tauOld
  iterVE <- 0;  stopVE <- 0

  while ((iterVE < currentOptions$maxIterVE) & (stopVE == 0)) {

    iterVE <- iterVE + 1


    tauOld <- tau;
    tau <- t(sapply(1:n, function(i){ # i <- 3
      indexListIFirst <- which(indexList[, 1] == i)
      indexListISecond <- which(indexList[, 2] == i)
      sapply(1:nbBlocks, function(k){ # k <- 1
        log(theta$blockProp[k]) +
          sum(logA[indexListIFirst, k, ] * tauOld[indexList[indexListIFirst, 2], ]) +
          sum(logA[indexListISecond, , k] * tauOld[indexList[indexListISecond, 1], ])
      })
    }))
    tau <- tau - apply(tau, 1, max)
    tau <- exp(tau); tau <- tau / rowSums(tau)
    tau <- tau + currentOptions$tauTol; tau <- tau / rowSums(tau)
    dTau <- distTau(tau,tauOld)
    if (dTau < currentOptions$valStopCrit)   {stopVE <- 1}
    #print(c(iterVE,dTau))
    if (iterVE == currentOptions$maxIterVE) {noConvergence = noConvergence + 1}
  }

  # tau

  # psi
  psi <- matrix(0, N, 2)
  sapply(1:nbBlocks, function(k){sapply(1:nbBlocks, function(l){ # k <- 1; l <- 2
    psi[, 2] <<- psi[, 2] + eta[, k, l] * tau[indexList[, 1], k] *  tau[indexList[, 2], l]
  })})
  psi[, 1]  = 1 - psi[, 2]

  qDist <- list(eta = eta, tau = tau,psi = psi)
  return(qDist)
}

#-------------------------------------------------------------------------------
# Computation of the lowerbound
#-------------------------------------------------------------------------------
lowerBoundScoreSBM <- function(scoreMat,theta,qDist,directed){

  # scoreMat <- scoreMat; theta <- thetaHat; qDist <- qDist

  # Dimensions
  nbBlocks <- length(theta$blockProp);
  N <- nrow(qDist$eta); n <- nbPairs2n(N, symmetric=!directed)
  indexList <- indices(n, symmetric=!directed)

  # Blocks
  espLogpZ <- sum(qDist$tau%*%log(theta$blockProp))
  HqZ <- -sum(qDist$tau*log(qDist$tau + (qDist$tau==0)))

  # Network
  tauArray <- array(dim=c(N, nbBlocks, nbBlocks))
  sapply(1:nbBlocks, function(k){sapply(1:nbBlocks, function(l){
    tauArray[, k, l] <<- qDist$tau[indexList[, 1], k] * qDist$tau[indexList[, 2], l]
  })})
  logConnectParam <- log(theta$connectParam); log1_ConnectParam <- log(1-theta$connectParam)
  espLogpG <- sum(sapply(1:nbBlocks, function(k){sapply(1:nbBlocks, function(l){
    sum(tauArray[, k, l] * (
      qDist$eta[, k, l] * logConnectParam[k, l] + (1 - qDist$eta[, k, l]) * log1_ConnectParam[k, l]
    ))
  })}))
  HqG <- sum(sapply(1:nbBlocks, function(k){sapply(1:nbBlocks, function(l){
    sum(tauArray[, k, l] * (
      qDist$eta[, k, l] * log(qDist$eta[, k, l] + (qDist$eta[, k, l]==0)) +
        (1 - qDist$eta[, k, l]) * log(1 - qDist$eta[, k, l] + (qDist$eta[, k, l]==1))
    ))
  })}))

  # Scores
  espLogpS <- sum(sapply(1:nbBlocks, function(k){sapply(1:nbBlocks, function(l){
    sum(tauArray[, k, l] * (
      (1 - qDist$eta[, k, l]) * theta$logPhi[, 1] + qDist$eta[, k, l] * theta$logPhi[, 2]
    ))
  })}))

  # Entropy and lower bound
  entropy <- HqZ + HqG
  lowerBound <- espLogpZ + espLogpG + espLogpS + entropy
  res <- list(espLogpZ=espLogpZ, HqZ=HqZ, espLogpG=espLogpG, HqG=HqG, espLogpS=espLogpZ,
              klZ=-espLogpZ-HqZ, klG=-espLogpG-HqG,
              entropy=entropy, lowerBound=lowerBound)
  return(res)

}




#-------------------------------------------------------------------------------
#----------  VEM algorithm ----------------------------------------------------
#-------------------------------------------------------------------------------

VEMScoreSBM <- function(scoreMat, directed, init,estimOptions = list(),monitoring = list(lowerBound = FALSE)){


  currentOptions <- list(
    verbosity = 0,
    maxIterVE = 100 ,
    maxIterVEM = 1000,
    tauTol = 2 * .Machine$double.eps,
    valCritStop = 1e-6,
    etaTol = 2 * .Machine$double.eps
  )
  currentOptions[names(estimOptions)] <- estimOptions



  #------------------------------------------------------------

  iterVEM <- 0
  deltaTau <- Inf

  if (monitoring$lowerBound) J  <- numeric(currentOptions$maxIterVEM)
  qDist = init



  #--------------   Algo begins


  while ((iterVEM < currentOptions$maxIterVEM) & (deltaTau > currentOptions$valCritStop))
  {

    iterVEM <- iterVEM + 1
    tauCurrent <- qDist$tau

    if (currentOptions$verbosity > 0){
      print(c(iterVEM,deltaTau))
    }
    #   M step ------------------
    theta <- mStepScoreSBM(scoreMat, qDist, directed)


    #  VE step ----------------
    qDist <- veStepScoreSBM(scoreMat, theta,tauOld = qDist$tau, directed,currentOptions)
    if (monitoring$lowerBound) J[iterVEM] =  lowerBoundScoreSBM(scoreMat,theta,qDist,directed)$lowerBound

    #  Stop check  ----------------
    deltaTau <- distTau(tauCurrent, qDist$tau)


  }
  #--------------  End of the algorithm
  #-------------  Reorder

  ord <- order(diag(theta$connectParam), decreasing  = TRUE)
  theta$blockProp <- theta$blockProp[ord]
  theta$connectParam <- theta$connectParam[ord,ord]

  output <- list(tau  = qDist$tau[,ord], theta = theta)
  if (monitoring$lowerBound) output$lowerBound <- J[1:iterVEM]

  return(output)
}


#-------------------------------------------------------------
#' Inference of the Score SBM  Model
#------------------------------------------------------------
#' \code{estimateScoreSBM} performs the estimation and model selection for the ScoreSBM model.
#' @param scoreList   :  list of Scores for each dyad  of an underlying network
#' @param directed    :  if true the inference network is directed. Default value = FALSE.
#' @param estimOptions : tunes the optimization process (see details below)
#' @param monitoring : specifies if the lowerBound along the VEM iterations is saved (monitoring = list(lowerBound = TRUE))
#' @details  The list of parameters \code{estimOptions} essentially tunes the optimization process and the variational EM algorithm, with the following parameters
#'  \itemize{
#'  \item{"verbosity"}{ controls the verbosity of the procedure (0 or 1). Default is 1.}
#'  \item{"exploreFactor"}{ controls the exploration of the number of groups in the initialization step, Default is 1.5}
#'  \item{"nbBlocksRange"}{ minimal and maximal number or blocks explored. Default is c(1,Inf)}
#'  \item{"nbCores"}{ integer for number of cores used. Default is 1. }
#'  \item{"maxIterVE"}{ Maximum number of iterations in the VE step. Default is 100.}
#'  \item{"maxIterVE"}{ Maximum number of iterations in the VEM algorithm. Default is 1000.}
#'  \item{"tauTol"}{ Tolerance in the VEM algorithm. Default value \code{tauTol = 2 * .Machine$double.eps}}
#'  \item{"etaTol"}{ Tolerance in the VEM algorithm. Default value \code{etaTol = 2 * .Machine$double.eps}}
#'  \item{"valCritStop"}{ Algorithm stops when the difference on successive tau is below valCritStop. Default value is  1e-6}
#'  }

#' @return The output is a list of estimated models, each one corresponding to a number of blocks.
#' @details Each element of the output list contains the following quantites:
#'  \itemize{
#'  \item{"theta"}{ estimated parameters}
#'  \item{"nbBlocks"}{ number of blocks in the underlying network}
#'  \item{"qDist"}{ variational approximation of the the conditional distribution q(G,Z | data)}
#'  \item{"ICL": }{Integrated Likelihood Criterion}
#'  \item{"pen": }{Penalty term for model selection (included in the ICL)}
#'  \item{"lowerBound": }{Lowerbound along the VEM iterations (provided if \code{monitoring$lowerBound = TRUE})}
#'  }
#' @examples
#' nbNodes  <- 60
#' directed <- FALSE
#' blockProp <- c(1/3,1/2,1/6)
#' nbBlocks   <- length(blockProp)
#' connectParam <- matrix(rbeta(nbBlocks^2,1.5,1.5 ),nbBlocks,nbBlocks)
#' connectParam <- 0.5*(connectParam + t(connectParam))
#' emissionParam <- list()
#' nbScores <- 4
#' emissionParam$noEdgeParam <- list(mean=rep(0,nbScores));
#' emissionParam$noEdgeParam$var <- diag(0.1,nrow = nbScores,ncol = nbScores)
#' emissionParam$EdgeParam <- list( mean= 1:nbScores)
#' emissionParam$EdgeParam$var <-  diag(0.1,nrow = nbScores,ncol = nbScores)
#' dataSim <- rScoreSBM(nbNodes,directed = TRUE, blockProp,connectParam,emissionParam,seed = NULL)
#' scoreList <- dataSim$ScoreNetworks
#' resEstim <- estimateScoreSBM(scoreList,directed)
#' @importFrom parallel mclapply
#' @importFrom pbmcapply pbmclapply



optimizeScoreSBM = function(scores,directed = FALSE, estimOptions=list(), monitoring = list()){

  currentOptions <- list(
    verbosity     = 1,
    explorFactor  = 1.5,
    nbBlocksRange = c(1,Inf),
    nbCores       = 1,
    maxIterVE = 100 ,
    maxIterVEM = 1000,
    tauTol = 2 * .Machine$double.eps,
    valCritStop = 1e-6,
    etaTol = 2 * .Machine$double.eps
  )
  currentMonitoring = list(
    lowerBound = TRUE
  )

  ## Current options are default expect for those passed by the user
  currentOptions[names(estimOptions)] <- estimOptions
  currentMonitoring[names(monitoring)] <- monitoring


  #------------------------ transform Data
  nbScores <- length(scores)
  scoreMat <- sapply(1:nbScores , function(q) {mat2Vect(scores[[q]], symmetric = !directed, diag = F)})

  #--------------------------------- initialisation .
  if (currentOptions$verbosity > 0) { print("-------------- Initialization ----------- ")}
  initAll <- initInferenceScoreSBM(scores, directed,currentOptions)

  #-------------------------------------------------------------------------
  #--------------------ESTIMATION ALL MODELS -------------------------------
  #------------------------------------------------------------------------
  nbModels <- length(initAll$tau)
  nbBlocksAll = vapply(1:nbModels,function(m){ncol(initAll$tau[[m]])},1)
  indModels <- (1:nbModels)[nbBlocksAll <= currentOptions$nbBlocksRange[2]]
  indModels <- (1:nbModels)[nbBlocksAll >= currentOptions$nbBlocksRange[1]]
  nbBlocksAll <- nbBlocksAll[indModels]
  if (currentOptions$verbosity > 0) { print("-------------- Estimation of the models ----------- ")}

  init_m  <- list(psi = initAll$psi)
  Estim_m <- function(m){
    if ((currentOptions$verbosity > 0 ) & (currentOptions$nbCores == 1)) {
      print(paste("-------------- Estimation for",m ,"blocks---- ----- ",sep = ' '))
    }
    init_m$tau = initAll$tau[[m]]
    init_m$eta = initAll$eta[[m]]
    resVEM_m <- VEMScoreSBM(scoreMat, directed, init_m,currentOptions,currentMonitoring)
    if ((currentOptions$verbosity > 0 ) & (currentOptions$nbCores == 1)) {
      print(paste("ICL = ",resVEM_m$ICL,sep = ' '))
    }

    return(resVEM_m)
  }

  listRes <- pbmclapply(indModels,Estim_m,mc.cores = currentOptions$nbCores);

  #-------------------------------------------------------------------------
  #--------------------REORDERING  MODELS -------------------------------
  #------------------------------------------------------------------------
  ordModels  <- order(sapply(indModels,function(m){listRes[[m]]$ICL}),decreasing  = TRUE)
  res <- lapply(ordModels, function(m){listRes[[m]]})
  return(res)

}

