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

test_that("Construction, fields access and other basics work in class ScoreSBM_Sampler (O/1 network observed through scores)", {

  ## SIMPLE UNDIRECTED BERNOULLI SBM with 3 scores


  ## Basic construction - check for wrong specifications
  mySampler <- ScoreSBM_sampler$new(nbNodes, FALSE, blockProp, connectParam, emissionParam)


  emissionParam2 <- emissionParam;
  emissionParam2$noEdgeParam$sigma2 <- -emissionParam$noEdgeParam$sigma2
  expect_error(ScoreSBM_sampler$new(nbNodes, FALSE, blockProp, connectParam, emissionParam2))


  emissionParam3 <- emissionParam;
  emissionParam3$noEdgeParam$mu <- c(-1,0)
  expect_error(ScoreSBM_sampler$new(nbNodes, FALSE, blockProp, connectParam,  emissionParam3))


  emissionParam4 <- emissionParam;
  emissionParam4$edgeParam$sigma2 <- matrix(0.1, nbScores+1,nbScores)
  expect_error(ScoreSBM_sampler$new(nbNodes, FALSE, blockProp, connectParam,  emissionParam3))


  ## Checking class
  expect_true(inherits(mySampler, "SBM"))
  expect_true(inherits(mySampler, "SBM_sampler"))
  expect_true(inherits(mySampler, "SimpleSBM_sampler"))
  expect_true(inherits(mySampler, "ScoreSBM_sampler"))
  ## Checking field access and format

  ## parameters
  expect_true(all(is.na(diag(mySampler$netMatrix))))
  expect_false(mySampler$directed)

  expect_true(isSymmetric(mySampler$scores[[1]]))
  expect_true(isSymmetric(mySampler$scores[[2]]))
  expect_true(isSymmetric(mySampler$scores[[3]]))

  ## blocks




})


