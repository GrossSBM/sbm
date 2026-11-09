# Estimation for multipartite SBM

Estimation for multipartite SBM

## Usage

``` r
estimateMultipartiteSBM(listSBM, estimOptions = list())
```

## Arguments

- listSBM:

  list of networks that were defined by the `defineSBM` function

- estimOptions:

  options for the inference procedure

## Value

a MultipartiteSBM_fit object with the estimated parameters and the
blocks in each Functional Group

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

## Missing values

The networks may contain missing values (`NA`): the corresponding dyads
are considered as not observed and are ignored in the inference
(performed by the package GREMLINS). This assumes that the dyads are
missing at random.

## Examples

``` r
if (FALSE) { # \dontrun{
## About the Parts/Functional Groups (FG)
blockProp <- list(c(0.16 ,0.40 ,0.44),c(0.3,0.7)) # prop of blocks in each FG
archiMultipartite <-  rbind(c(1,2),c(2,2),c(1,1)) # architecture of the multipartite net.
nbNodes <- c(60,50)
## About the connection matrices
directed <- c(NA, TRUE, FALSE) # type of each network
model <- c('gaussian','bernoulli','poisson')
C1 <-
 list(mean = matrix(c(6.1, 8.9, 6.6, 9.8, 2.6, 1.0), 3, 2),
      var  = matrix(c(1.6, 1.6, 1.8, 1.7 ,2.3, 1.5),3, 2))
C2 <- list(mean = matrix(c(0.7,1.0, 0.4, 0.6),2, 2))
m3 <- matrix(c(2.5, 2.6 ,2.2 ,2.2, 2.7 ,3.0 ,3.6, 3.5, 3.3),3,3 )
C3 <- list(mean = .5 * (m3 + t(m3)))
connectParam <- list(C1, C2, C3)
## Graph Sampling
mySampleMSBM <- sampleMultipartiteSBM(nbNodes, blockProp,
                                      archiMultipartite, connectParam, model,
                                      directed, dimLabels = c('A','B'), seed = 2)
listSBM <- mySampleMSBM$listSBM
estimOptions <- list(initBM = FALSE, nbCores  = 2)
myMSBM <- estimateMultipartiteSBM(listSBM, estimOptions)
plot(myMSBM, type = "data")
plot(myMSBM, type = "expected")
plot(myMSBM, type = "meso")
} # }
```
