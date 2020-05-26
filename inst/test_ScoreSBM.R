## Common parameters

nbNodes  <- 60
directed <- FALSE
blockProp <- c(1/3,1/2,1/6)
nbBlocks   <- length(blockProp)
connectParam <- matrix(rbeta(nbBlocks^2,1.5,1.5 ),nbBlocks,nbBlocks)
connectParam <- list(mu = 0.5*(connectParam + t(connectParam)))
emissionParam <- list()
nbScores <- 4
emissionParam$noEdgeParam <- list(mu=rep(0,nbScores));
emissionParam$noEdgeParam$sigma2 <- diag(0.1,nrow = nbScores,ncol = nbScores)
emissionParam$edgeParam <- list( mu= 1:nbScores)
emissionParam$edgeParam$sigma2 <-  diag(0.1,nrow = nbScores,ncol = nbScores)

mySampler <- ScoreSBM_sampler$new(nbNodes, FALSE, blockProp, connectParam, emissionParam,seed = 3)

S <- sampleScoreSBM(nbNodes,blockProp,connectParam, directed = FALSE,emissionParam, seed = NULL)
scoresList <- S$scores
length(scoresList)


estim <- ScoreSBM_fit$new(scoresList, directed = FALSE)
length(estim$scores)

estim$optimize()

library(NoisySBM)
res1 <- NoisySBM::estimateNoisySBM(scoresList,directed = FALSE, estimOptions=list(), monitoring = list())
