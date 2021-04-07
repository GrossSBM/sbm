Nnodes <- c(40,30)
blockProp <- list(c(.4,.6),c(0.5,0.5))
nbLayers <- 3
connectParam <- list(list(mean=matrix(rbeta(4,.5,.5),2,2)),list(mean=matrix(rexp(4,.5),2,2)))
connectParam[[3]] <- list(mean=matrix(rbeta(4,0.9,.9),2,2))
model <- c("bernoulli","poisson","bernoulli")
type <- "bipartite"
mySampleMultiplexSBM <-
  SampleMultiplexSBM(
    nbNodes = Nnodes,
    blockProp = blockProp,
    nbLayers = nbLayers,
    connectParam = connectParam,
    model=model,
    type=type)
listSBM <- mySampleMultiplexSBM$listSBM
fitMultiplex <- estimateMultiplexSBM(listSBM)


MultipartiteSBM_fit$new(listSBM)
plotMyMultipartiteMatrix(listSBM)



listSBM2 <- listSBM
listSBM2[[1]]$dimLabels = 'robert'
listSBM2[[1]]$dimLabels
myMSBMObject2 <- MultipartiteSBM_fit$new(listSBM2)
myMSBMObject2$architecture
myMSBMObject2$nbNodes
listSBM2[[1]]$nbNodes
myMSBMObject2$dimLabels
myMSBMObject2$nbNodes



listSBM3 <- listSBM
listSBM3[[1]]$dimLabels = c('actor','actorLayer1')
listSBM2[[1]]$dimLabels
myMSBMObject2 <- MultipartiteSBM_fit$new(listSBM2)
myMSBMObject2$architecture
myMSBMObject2$nbNodes
listSBM2[[1]]$nbNodes
myMSBMObject2$dimLabels
myMSBMObject2$nbNodes

