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
estimOptions <- list(initBM = FALSE, nbCores  = 2)

myMultiplexSBM <- estimateMultiplexSBM(listSBM,estimOptions,dependent=TRUE)
myMultiplexSBM$connectParam



Q <- 3
nbLayers <- 2
connectParam <- list()
connectParam <- vector("list",nbLayers)
connectParam[[1]]$mean <-  matrix(.1,Q,Q) + diag(1:Q)
connectParam[[2]]$mean <- matrix(-2,Q,Q) + diag(rev(Q:1))
connectParam[[1]]$var <- matrix(c(2,1,1,4),Q,Q)
connectParam[[2]]$var <- matrix(c(2,1,1,4),Q,Q)

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
    dependent=F)
listSBM <- mySampleMultiplexSBM$listSBM
estimOptions <- list(initBM = FALSE, nbCores  = 2)

myMultiplexSBM <- estimateMultiplexSBM(listSBM,estimOptions,dependent=F)
myMultiplexSBM$connectParam


