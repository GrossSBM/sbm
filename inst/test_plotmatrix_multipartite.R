
nbBlocks  <- c(3,2,2)
nbNodes <-  c(100,50,40)
blockProp <- vector("list", 3)
blockProp[[1]] <- c(0.4,0.3,0.3) # in Functional Group 1
blockProp[[2]] <- c(0.6,0.4) # in Functional Group 2
blockProp[[3]]  <- c(0.6,0.4) # in Functional Group 3
archiMultipartite  <-  rbind(c(1,2),c(2,3),c(2,2),c(1,3))
model <- c('bernoulli','poisson','bernoulli','gaussian')
directed <- c( NA, NA  ,  FALSE , NA)
connectParam <- list()
E <- archiMultipartite
mu <- rbeta(nbBlocks[E[1,1]] * nbBlocks[E[1,2]],1,1 )
connectParam[[1]] <- list(mean = matrix(mu,nrow = nbBlocks[E[1,1]], ncol = nbBlocks[E[1,2]] ))
mu <- rgamma(nbBlocks[E[2,1]] * nbBlocks[E[2,2]],7.5,0.01 )
connectParam[[2]] <- list(mean  =  matrix(mu,nrow = nbBlocks[E[2,1]], ncol = nbBlocks[E[2,2]]))
p <- rbeta(nbBlocks[E[3,1]] * nbBlocks[E[3,2]],0.9,0.9 )
p <- 1/2*(p + t(p))
connectParam[[3]] <- list(mean  =  matrix(p, nrow = nbBlocks[E[3,1]], ncol = nbBlocks[E[3,2]]))
mu <- rnorm(nbBlocks[E[4,1]] * nbBlocks[E[4,2]],7.5,10 )
connectParam[[4]] <- list(mean = matrix(mu, nrow = nbBlocks[E[4,1]], ncol = nbBlocks[E[4,2]]))
v <- rgamma(nbBlocks[E[4,1]] * nbBlocks[E[4,2]],7.5,0.1 )
connectParam[[4]]$var <- matrix(v, nrow = nbBlocks[E[4,1]], ncol = nbBlocks[E[4,2]])
## Graph Sampling
mySampleMSBM <- sampleMultipartiteSBM(nbNodes, blockProp,
                                      archiMultipartite, connectParam, model,
                                      directed, dimLabels = as.list(c('A','B','C')))
listSBM <- mySampleMSBM$listSBM
estimOptions = list(initBM = FALSE)
myMSBM <- estimateMultipartiteSBM(listSBM,estimOptions)

myMSBM
listMat <- lapply(myMSBM$listSBM,function(s){s$netMatrix})
E <- myMSBM$archiMultipartite
namesFG <- myMSBM$dimLabels
nbNodes <- myMSBM$nbNodes
normalized <- TRUE
clustering <- myMSBM$memberships
distrib <- myMSBM$modelName


plotMultipartiteMatrix(listMat, E, nbNodes, namesFG, normalized, distrib, clustering)
