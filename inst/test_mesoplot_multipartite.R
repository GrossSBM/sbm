myMSBM <- myEstimatedMSBM
E <- myMSBM$archiMultipartite
theta <- myMSBM$connectParam
list_pi <- myMSBM$blockProp
v_distrib <- myMSBM$modelName
nbNodes <- myMSBM$nbNodes
nodeLabels <- myMSBM$dimLabels
directed<- myMSBM$directed
plotOptions <- NULL

P <- plotMesoMultipartite(E,theta, list_pi,v_distrib,directed,nbNodes,nodeLabels,plotOptions)




nbBlocks  <- c(3,2,2) #number of clusters in each functional group
nbNodes <-  c(100,50,40)
blockProp <- vector("list", 3)  # parameters of clustering in each functional group
blockProp[[1]] <- c(0.4,0.3,0.3) # in Functional Group 1
blockProp[[2]] <- c(0.6,0.4) # in Functional Group 2
blockProp[[3]]  <- c(0.6,0.4) # in Functional Group 3
# About the interactions between the FG
archiMultipartite  <-  rbind(c(1,2),c(2,3),c(2,2),c(1,3)) #
model <- c('bernoulli','poisson','bernoulli','gaussian') # type of distribution in each network
# for each network : directed or not (not required for an interaction between two different FG)
directed <- c( NA, NA  ,  TRUE , NA)
connectParam <- list()
E <- archiMultipartite
m <- rbeta(nbBlocks[E[1,1]] * nbBlocks[E[1,2]],1,1 )
connectParam[[1]] <- list(mean = matrix(m,nrow = nbBlocks[E[1,1]], ncol = nbBlocks[E[1,2]] ))
m <- rgamma(nbBlocks[E[2,1]] * nbBlocks[E[2,2]],7.5,0.01 )
connectParam[[2]] <- list(mean  =  matrix(m,nrow = nbBlocks[E[2,1]], ncol = nbBlocks[E[2,2]]))
p <- rbeta(nbBlocks[E[3,1]] * nbBlocks[E[3,2]],0.9,0.9 )
#p <- 1/2*(p + t(p))
connectParam[[3]] <- list(mean  =  matrix(p, nrow = nbBlocks[E[3,1]], ncol = nbBlocks[E[3,2]]))
m <- rnorm(nbBlocks[E[4,1]] * nbBlocks[E[4,2]],7.5,10 )
connectParam[[4]] <- list(mean = matrix(m, nrow = nbBlocks[E[4,1]], ncol = nbBlocks[E[4,2]]))
v <- rgamma(nbBlocks[E[4,1]] * nbBlocks[E[4,2]],7.5,0.1 )
connectParam[[4]]$var <- matrix(v, nrow = nbBlocks[E[4,1]], ncol = nbBlocks[E[4,2]])
dimLabels <- as.list(c('A','B','C'))
seed <- 10
## Graph Sampling
mySampleMBM <- sampleMultipartiteSBM(nbNodes, blockProp, archiMultipartite,
                                     connectParam, model, directed, dimLabels,seed)

myEstimatedMSBM <- estimateMultipartiteSBM(mySampleMBM$listSBM)
plot(myEstimatedMSBM,normalized = TRUE)
plot(myEstimatedMSBM,type='meso')
