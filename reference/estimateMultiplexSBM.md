# Estimation for Multiplex SBM

Estimation for Multiplex SBM

## Usage

``` r
estimateMultiplexSBM(listSBM, dependent = FALSE, estimOptions = list())
```

## Arguments

- listSBM:

  list of networks that were defined by the `defineSBM` function

- dependent:

  logical parameter indicating whether the networks in the multiplex
  structure are dependent beyond the latent variables,

- estimOptions:

  options for the inference procedure

## Value

a MultiplexSBM_fit object with the estimated parameters and the blocks

## Details

The list of parameters `estimOptions` essentially tunes the optimization
process and the variational EM algorithm, with the following parameters

- "nbCores" integer for number of cores used. Default is 2 (1 under
  Positron, which does not allow forking)

- "verbosity" integer for verbosity (0, 1). Default is 1

- "nbBlocksRange" List of length the number of functional groups, each
  element supplying the minimal and maximal number of blocks to be
  explored. The names of the list must be the names of the functional
  groups. Default value is from 1 to 10)

- "initBM" Boolean. True if using simple and bipartite SBM as
  initialisations. Default value = TRUE

- "maxiterVEM" Number of max. number of iterations in the VEM. Default
  value = 100

- "maxiterVE" Number of max. number of iterations in the VE. Default
  value = 100

- "plot" boolean, should the ICL be dynamically plotted or not. Default
  is TRUE, FALSE when "verbosity" is 0. For dependent networks

- "exploreFactor" control the exploration of the number of groups. For
  dependent networks

- "exploreMin" explore at least until exploreMin even if the exploration
  factor rule is achieved. Default 4. See the package blockmodels for
  details. For dependent networks

- "exploreMax" Stop exploration at exploreMax even if the exploration
  factor rule is not achieved. Default Inf. See the package blockmodels
  for details. For dependent networks

- "nbBlocksRange" minimal and maximal number or blocks explored. For
  dependent networks

- "fast" logical: should approximation be used for Bernoulli model with
  covariates. Default to `TRUE`. For dependent networks

## Examples

``` r
if (FALSE) { # \dontrun{
### =======================================
### MULTIPLEX SBM without dependence between layers
##
Nnodes <- 40
blockProp <- c(.4,.6)
nbLayers <- 2
connectParam <- list(list(mean=matrix(rbeta(4,.5,.5),2,2)),list(mean=matrix(rexp(4,.5),2,2)))
model <- c("bernoulli","poisson")
type <- "directed"
mySampleMultiplexSBM <-
   sampleMultiplexSBM(
   nbNodes = Nnodes,
    blockProp = blockProp,
   nbLayers = nbLayers,
   connectParam = connectParam,
   model=model,
   type=type)
listSBM <- mySampleMultiplexSBM$listSBM
estimOptions <- list(initBM = FALSE, nbCores  = 2)
myMultiplexSBM <- estimateMultiplexSBM(listSBM,estimOptions,dependent=FALSE)
### =======================================
### MULTIPLEX SBM Gaussian with dependence
##
Q <- 3
nbLayers <- 2
connectParam <- list()
connectParam$mu <- vector("list",nbLayers)
connectParam$mu[[1]] <-  matrix(.1,Q,Q) + diag(1:Q)
connectParam$mu[[2]] <- matrix(-2,Q,Q) + diag(rev(Q:1))
connectParam$Sigma <- matrix(c(2,1,1,4),nbLayers,nbLayers)
model <- rep("gaussian",2)
type <- "directed"
Nnodes <- 80
blockProp <- c(.3,.3,.4)
mySampleMultiplexSBM <-
  sampleMultiplexSBM(
     nbNodes = Nnodes,
     blockProp = blockProp,
     nbLayers = nbLayers,
     connectParam = connectParam,
     model=model,
     type="undirected",
     dependent=TRUE)
listSBM <- mySampleMultiplexSBM$listSBM
myMultiplexSBM <- estimateMultiplexSBM(listSBM,estimOptions,dependent=TRUE)
## MultiplexSBM Bernoulli with dependence
Q <- 2
P00<-matrix(runif(Q*Q),Q,Q)
P10<-matrix(runif(Q*Q),Q,Q)
P01<-matrix(runif(Q*Q),Q,Q)
P11<-matrix(runif(Q*Q),Q,Q)
SumP<-P00+P10+P01+P11
P00<-P00/SumP
P01<-P01/SumP
P10<-P10/SumP
P11<-P11/SumP
connectParam <- list()
connectParam$prob00 <- P00
connectParam$prob01 <- P01
connectParam$prob10 <- P10
connectParam$prob11 <- P11
model <- rep("bernoulli",2)
type <- "directed"
nbLayers <- 2
Nnodes <- 40
blockProp <- c(.6,.4)
mySampleMultiplexSBM <-
   sampleMultiplexSBM(
     nbNodes = Nnodes,
     blockProp = blockProp,
     nbLayers = nbLayers,
     connectParam = connectParam,
     model=model,
     type=type,
     dependent=TRUE)
listSBM <- mySampleMultiplexSBM$listSBM
myMultiplexSBM <- estimateMultiplexSBM(listSBM,estimOptions,dependent=TRUE)
} # }
```
