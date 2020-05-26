## Common parameters
nbNodes  <- 90
nbBlocks <- 3
blockProp <- c(.5, .25, .25) # group proportions
means <- diag(.4, 3) + 0.05
connectParam <- list(mu = means)
nbScores <- 3;
emissionParam <- list(edgeParam = list(),noEdgeParam = list())
emissionParam$noEdgeParam$mu <- c(-1,0,1)
emissionParam$edgeParam$mu <- c(10,13,15)
emissionParam$noEdgeParam$sigma2 <- matrix(0.1, nbScores,nbScores);
diag(emissionParam$noEdgeParam$sigma2)  = 1;
emissionParam$edgeParam$sigma2 <- matrix(1, nbScores,nbScores)
diag(emissionParam$edgeParam$sigma2)  = 2;

mySampler <- ScoreSBM_sampler$new(nbNodes, FALSE, blockProp, connectParam, emissionParam,seed = 3)

S <- sampleScoreSBM(nbNodes,blockProp,connectParam, directed = FALSE,emissionParam, seed = NULL)
scoresList <- S$scores
length(scoresList)


estim <- ScoreSBM_fit$new(scoresList, directed = FALSE)
length(estim$scores)

estim$optimize()
