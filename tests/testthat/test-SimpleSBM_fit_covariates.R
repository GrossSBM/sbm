set.seed(1234)

rmse <- function(theta, theta_star) { sqrt(sum((theta - theta_star)^2)/sum(theta_star^2)) }

## Common parameters
nbNodes  <- 40
nbBlocks <- 2
blockProp <- c(.5, .5) # group proportions
covarParam <- c(-2,2)
dimLabels <- list(row = "rowLabel", col = "colLabel")
covar1 <- matrix(rnorm(nbNodes**2), nbNodes, nbNodes)
covar2 <- matrix(rnorm(nbNodes**2), nbNodes, nbNodes)
covarList_directed <- list(covar1 = covar1, covar2 = covar2)

covar1 <- covar1 + t(covar1)
covar2 <- covar2 + t(covar2)
covarList <- list(covar1 = covar1, covar2 = covar2)

test_that("SimpleSBM_fit 'Bernoulli' model, undirected, one covariate", {

  ## SIMPLE UNDIRECTED BERNOULLI SBM
  means <- diag(.4, 2) + 0.05
  connectParam <- list(mean = means)

  ## Basic construction - check for wrong specifications
  mySampler <- SimpleSBM$new('bernoulli', nbNodes, FALSE, blockProp, connectParam, covarParam = covarParam[1], covarList = covarList[1])
  mySampler$rMemberships(store = TRUE)
  mySampler$rEdges(store = TRUE)

  ## Construction----------------------------------------------------------------
  mySBM <- SimpleSBM_fit$new(mySampler$networkData, 'bernoulli', FALSE, covarList = covarList[1])
  expect_error(SimpleSBM_fit$new(mySampler$networkData, 'bernouilli', FALSE, covarList = covarList[1]))
  expect_error(SimpleSBM_fit$new(mySampler$networkData[1:20, 1:30], 'bernouilli', FALSE, covarList = covarList[1]))
  expect_error(SimpleSBM_fit$new(mySampler$networkData, 'bernoulli', TRUE, covarList = covarList[1]))
  expect_error(SimpleSBM_fit$new(mySampler$networkData, 'bernoulli', FALSE, covarList = covarList[[1]]))

  ## Checking class
  expect_true(inherits(mySBM, "SBM"))
  expect_true(inherits(mySBM, "SimpleSBM"))
  expect_true(inherits(mySBM, "SimpleSBM_fit"))

  ## Checking field access and format prior to estimation
  ## parameters
  expect_equal(mySBM$modelName, 'bernoulli')
  expect_equal(unname(mySBM$nbNodes), nbNodes)
  expect_equal(mySBM$dimLabels, c(node="nodeName"))
  expect_equal(mySBM$nbDyads, nbNodes*(nbNodes - 1)/2)
  expect_true(all(is.na(diag(mySBM$networkData))))
  expect_true(isSymmetric(mySBM$networkData))
  expect_true(!mySBM$directed)
  expect_true(is.matrix(mySBM$connectParam$mean))

  ## covariates
  expect_true(all(dim(mySBM$covarEffect) == c(nbNodes, nbNodes)))
  expect_equal(mySBM$nbCovariates, 1)
  expect_equal(mySBM$covarList, covarList[1])
  expect_equal(mySBM$covarParam, c(0))

  ## S3 methods
  expect_equal(coef(mySBM, 'connectivity'), mySBM$connectParam)
  expect_equal(coef(mySBM, 'block')       , mySBM$blockProp)
  expect_equal(coef(mySBM, 'covariates')  , mySBM$covarParam)

  ## Estimation-----------------------------------------------------------------
  BM_out <- mySBM$optimize(estimOptions = list(verbosity = 0, fast = TRUE))
  mySBM$setModel(2)

  ## Expectation
  expect_equal(dim(mySBM$expectation), c(nbNodes, nbNodes))
  expect_true(all(mySBM$expectation >= 0, na.rm = TRUE))
  expect_true(all(mySBM$expectation <= 1, na.rm = TRUE))
  expect_null(mySBM$connectParam$var)

  ## blocks
  expect_equal(mySBM$nbBlocks, nbBlocks)
  expect_equal(dim(mySBM$probMemberships), c(nbNodes, nbBlocks))
  expect_equal(sort(unique(mySBM$memberships)), 1:nbBlocks)
  expect_equal(length(mySBM$memberships), nbNodes)

  ## S3 methods
  expect_equal(coef(mySBM, 'connectivity'), mySBM$connectParam)
  expect_equal(coef(mySBM, 'block')       , mySBM$blockProp)
  expect_equal(coef(mySBM, 'covariates')  , mySBM$covarParam)
  expect_equal(mySBM$predict(), predict(mySBM))
  expect_equal(predict(mySBM), fitted(mySBM))
  expect_equal(predict(mySBM, covarList[1]), fitted(mySBM))
  expect_error(predict(mySBM, covarList))

  ## prediction wrt BM
  for (Q in mySBM$storedModels$indexModel) {
    pred_bm  <- BM_out$prediction(Q = Q)
    mySBM$setModel(Q)
    pred_sbm <- predict(mySBM)
    expect_lt( rmse(pred_bm, pred_sbm), 1e-12)
  }

})

test_that("SimpleSBM_fit 'Bernoulli' model, directed, one covariate", {

  ## SIMPLE UNDIRECTED BERNOULLI SBM
  means <- matrix(c(0.1, 0.4, 0.6, 0.9), 2,  2)
  connectParam <- list(mean = means)

  ## Basic construction - check for wrong specifications
  mySampler <- SimpleSBM$new('bernoulli', nbNodes, TRUE, blockProp, connectParam, c(node = "nodeName"), covarParam[1], covarList_directed[1])
  mySampler$rMemberships(store = TRUE)
  mySampler$rEdges(store = TRUE)

  ## Construction----------------------------------------------------------------
  mySBM <- SimpleSBM_fit$new(mySampler$networkData, 'bernoulli', TRUE, covarList = covarList_directed[1])
  expect_error(SimpleSBM_fit$new(mySampler$networkData, 'bernouilli', TRUE, covarList = covarList_directed[1]))
  expect_error(SimpleSBM_fit$new(mySampler$networkData, 'bernoulli', FALSE, covarList = covarList_directed[1]))
  expect_error(SimpleSBM_fit$new(mySampler$networkData[1:20, 1:30], 'bernouilli', FALSE, covarList = covarList_directed[1]))

  ## Checking class
  expect_true(inherits(mySBM, "SBM"))
  expect_true(inherits(mySBM, "SimpleSBM"))
  expect_true(inherits(mySBM, "SimpleSBM_fit"))

  ## Checking field access and format prior to estimation
  ## parameters
  expect_equal(mySBM$modelName, 'bernoulli')
  expect_equal(unname(mySBM$nbNodes), nbNodes)
  expect_equal(mySBM$dimLabels, c(node="nodeName"))
  expect_equal(mySBM$nbDyads, nbNodes*(nbNodes - 1))
  expect_true(all(is.na(diag(mySBM$networkData))))
  expect_true(!isSymmetric(mySBM$networkData))
  expect_true(mySBM$directed)
  expect_true(is.matrix(mySBM$connectParam$mean))

  ## covariates
  expect_true(all(dim(mySBM$covarEffect) == c(nbNodes, nbNodes)))
  expect_equal(mySBM$nbCovariates, 1)
  expect_equal(mySBM$covarList, covarList_directed[1])
  expect_equal(mySBM$covarParam, c(0))

  ## S3 methods
  expect_equal(coef(mySBM, 'connectivity'), mySBM$connectParam)
  expect_equal(coef(mySBM, 'block')       , mySBM$blockProp)
  expect_equal(coef(mySBM, 'covariates')  , mySBM$covarParam)

  ## Estimation-----------------------------------------------------------------
  BM_out <- mySBM$optimize(estimOptions  = list(verbosity = 0, fast = TRUE))
  mySBM$setModel(2)

  ## Expectation
  expect_equal(dim(mySBM$expectation), c(nbNodes, nbNodes))
  expect_true(all(mySBM$expectation >= 0, na.rm = TRUE))
  expect_true(all(mySBM$expectation <= 1, na.rm = TRUE))
  expect_null(mySBM$connectParam$var)

  ## blocks
  expect_equal(mySBM$nbBlocks, nbBlocks)
  expect_equal(dim(mySBM$probMemberships), c(nbNodes, nbBlocks))
  expect_equal(sort(unique(mySBM$memberships)), 1:nbBlocks)
  expect_equal(length(mySBM$memberships), nbNodes)

  ## S3 methods
  expect_equal(coef(mySBM, 'connectivity'), mySBM$connectParam)
  expect_equal(coef(mySBM, 'block')       , mySBM$blockProp)
  expect_equal(coef(mySBM, 'covariates')  , mySBM$covarParam)
  expect_equal(mySBM$predict(), predict(mySBM))
  expect_equal(predict(mySBM), fitted(mySBM))

  expect_equal(predict(mySBM, covarList[1]), fitted(mySBM))

  ## prediction wrt BM
  for (Q in mySBM$storedModels$indexModel) {
    pred_bm  <- BM_out$prediction(Q = Q)
    mySBM$setModel(Q)
    pred_sbm <- predict(mySBM)
    expect_lt( rmse(pred_bm, pred_sbm), 1e-12)
  }

})

test_that("SimpleSBM_fit 'Poisson' model, undirected, two covariates", {

  ## SIMPLE UNDIRECTED POISSON SBM
  means <- diag(15, 2) + 5
  connectParam <- list(mean = means)

  ## Basic construction - check for wrong specifications
  mySampler <- SimpleSBM$new('poisson', nbNodes, FALSE, blockProp, connectParam, covarParam = covarParam[1], covarList = covarList[1])
  mySampler$rMemberships(store = TRUE)
  mySampler$rEdges(store = TRUE)

  ## Construction----------------------------------------------------------------
  mySBM <- SimpleSBM_fit$new(mySampler$networkData, 'poisson', FALSE, covarList = covarList[1])
  expect_error(SimpleSBM_fit$new(SamplerBernoulli$networkData, 'poison', FALSE, covarList = covarList[1]))
  expect_error(SimpleSBM_fit$new(SamplerBernoulli$networkData[1:20, 1:30], 'poisson', FALSE, covarList = covarList[1]))
  expect_error(SimpleSBM_fit$new(SamplerBernoulli$networkData, 'poisson', TRUE, covarList = covarList[1]))

  ## Checking class
  expect_true(inherits(mySBM, "SBM"))
  expect_true(inherits(mySBM, "SimpleSBM"))
  expect_true(inherits(mySBM, "SimpleSBM_fit"))

  ## Checking field access and format prior to estimation
  ## parameters
  expect_equal(mySBM$modelName, 'poisson')
  expect_equal(unname(mySBM$nbNodes), nbNodes)
  expect_equal(mySBM$dimLabels, c(node="nodeName"))
  expect_equal(mySBM$nbDyads, nbNodes*(nbNodes - 1)/2)
  expect_true(all(is.na(diag(mySBM$networkData))))
  expect_true(isSymmetric(mySBM$networkData))
  expect_true(!mySBM$directed)
  expect_true(is.matrix(mySBM$connectParam$mean))

  ## covariates
  expect_true(all(dim(mySBM$covarEffect) == c(nbNodes, nbNodes)))
  expect_equal(mySBM$nbCovariates, 1)
  expect_equal(mySBM$covarList, covarList[1])
  expect_equal(mySBM$covarParam, c(0))

  ## S3 methods
  expect_equal(coef(mySBM, 'connectivity'), mySBM$connectParam)
  expect_equal(coef(mySBM, 'block')       , mySBM$blockProp)
  expect_equal(coef(mySBM, 'covariates')  , mySBM$covarParam)

  ## Estimation-----------------------------------------------------------------
  BM_out <- mySBM$optimize(estimOptions=list(verbosity = 0))
  mySBM$setModel(2)

  ## Expectation
  expect_equal(dim(mySBM$expectation), c(nbNodes, nbNodes))
  expect_true(all(mySampler$expectation >= 0, na.rm = TRUE))
  expect_null(mySBM$connectParam$var)

  ## blocks
  expect_equal(mySBM$nbBlocks, nbBlocks)
  expect_equal(dim(mySBM$probMemberships), c(nbNodes, nbBlocks))
  expect_equal(sort(unique(mySBM$memberships)), 1:nbBlocks)
  expect_equal(length(mySBM$memberships), nbNodes)

  ## S3 methods
  expect_equal(coef(mySBM, 'connectivity'), mySBM$connectParam)
  expect_equal(coef(mySBM, 'block')       , mySBM$blockProp)
  expect_equal(coef(mySBM, 'covariates')  , mySBM$covarParam)
  expect_equal(mySBM$predict(), predict(mySBM))
  expect_equal(predict(mySBM), fitted(mySBM))

  expect_equal(predict(mySBM, covarList[1]), fitted(mySBM))

  ## prediction wrt BM
  for (Q in mySBM$storedModels$indexModel) {
    pred_bm  <- BM_out$prediction(Q = Q)
    mySBM$setModel(Q)
    pred_sbm <- predict(mySBM)
    expect_lt( rmse(pred_bm, pred_sbm), 1e-12)
  }

})

test_that("SimpleSBM_fit 'Poisson' model, directed, two covariates", {

  ## SIMPLE UNDIRECTED BERNOULLI SBM
  means <- matrix(c(1, 4, 7, 9), 2,  2)
  connectParam <- list(mean = means)

  ## Basic construction - check for wrong specifications
  mySampler <- SimpleSBM$new('poisson', nbNodes, TRUE, blockProp, connectParam, covarParam = covarParam[1], covarList = covarList_directed[1])
  mySampler$rMemberships(store = TRUE)
  mySampler$rEdges(store = TRUE)

  ## Construction----------------------------------------------------------------
  mySBM <- SimpleSBM_fit$new(mySampler$networkData, 'poisson', TRUE, covarList = covarList_directed[1])
  expect_error(SimpleSBM_fit$new(SamplerBernoulli$networkData, 'poison', TRUE, covarList = covarList_directed[1]))
  expect_error(SimpleSBM_fit$new(SamplerBernoulli$networkData[1:20, 1:30], 'poisson', FALSE, covarList = covarList_directed[1]))
  expect_error(SimpleSBM_fit$new(SamplerBernoulli$networkData, 'poisson', FALSE, covarList = covarList_directed[1]))

  ## Checking class
  expect_true(inherits(mySBM, "SBM"))
  expect_true(inherits(mySBM, "SimpleSBM"))
  expect_true(inherits(mySBM, "SimpleSBM_fit"))

  ## Checking field access and format prior to estimation
  ## parameters
  expect_equal(mySBM$modelName, 'poisson')
  expect_equal(unname(mySBM$nbNodes), nbNodes)
  expect_equal(mySBM$dimLabels, c(node="nodeName"))
  expect_equal(mySBM$nbDyads, nbNodes*(nbNodes - 1))
  expect_true(all(is.na(diag(mySBM$networkData))))
  expect_true(!isSymmetric(mySBM$networkData))
  expect_true(mySBM$directed)
  expect_true(is.matrix(mySBM$connectParam$mean))

  ## covariates
    expect_true(all(dim(mySBM$covarEffect) == c(nbNodes, nbNodes)))
  expect_equal(mySBM$nbCovariates, 1)
  expect_equal(mySBM$covarList, covarList_directed[1])
  expect_equal(mySBM$covarParam, c(0))

  ## S3 methods
  expect_equal(coef(mySBM, 'connectivity'), mySBM$connectParam)
  expect_equal(coef(mySBM, 'block')       , mySBM$blockProp)
  expect_equal(coef(mySBM, 'covariates')  , mySBM$covarParam)

  ## Estimation-----------------------------------------------------------------
  BM_out <- mySBM$optimize(estimOptions=list(verbosity = 0))
  mySBM$setModel(2)

  ## Expectation
  expect_equal(dim(mySBM$expectation), c(nbNodes, nbNodes))
  expect_true(all(mySampler$expectation >= 0, na.rm = TRUE))
  expect_null(mySBM$connectParam$var)

  ## blocks
  expect_equal(mySBM$nbBlocks, nbBlocks)
  expect_equal(dim(mySBM$probMemberships), c(nbNodes, nbBlocks))
  expect_equal(sort(unique(mySBM$memberships)), 1:nbBlocks)
  expect_equal(length(mySBM$memberships), nbNodes)

  ## S3 methods
  expect_equal(coef(mySBM, 'connectivity'), mySBM$connectParam)
  expect_equal(coef(mySBM, 'block')       , mySBM$blockProp)
  expect_equal(coef(mySBM, 'covariates')  , mySBM$covarParam)
  expect_equal(mySBM$predict(), predict(mySBM))
  expect_equal(predict(mySBM), fitted(mySBM))

  expect_equal(predict(mySBM, covarList[1]), fitted(mySBM))

  ## prediction wrt BM
  for (Q in mySBM$storedModels$indexModel) {
    pred_bm  <- BM_out$prediction(Q = Q)
    mySBM$setModel(Q)
    pred_sbm <- predict(mySBM)
    expect_lt( rmse(pred_bm, pred_sbm), 1e-12)
  }

})


test_that("SimpleSBM_fit 'Gaussian' model, undirected, two covariates", {

  ## SIMPLE UNDIRECTED GAUSSIAN SBM
  means <- diag(15., 2) + 5 # connectivity matrix: affiliation network
  connectParam <- list(mean = means, var = 2)

  ## Basic construction - check for wrong specifications
  mySampler <- SimpleSBM$new('gaussian', nbNodes, FALSE, blockProp, connectParam, covarParam = covarParam, covarList = covarList)
  mySampler$rMemberships(store = TRUE)
  mySampler$rEdges(store = TRUE)

  ## Construction----------------------------------------------------------------
  mySBM <- SimpleSBM_fit$new(mySampler$networkData, 'gaussian', FALSE, covarList = covarList)
  expect_error(SimpleSBM_fit$new(SamplerBernoulli$networkData, 'normal', FALSE, covarList = covarList))
  expect_error(SimpleSBM_fit$new(SamplerBernoulli$networkData[1:20, 1:30], 'gaussian', FALSE, covarList = covarList))
  expect_error(SimpleSBM_fit$new(SamplerBernoulli$networkData, 'gaussian', TRUE, covarList = covarList))

  ## Checking class
  expect_true(inherits(mySBM, "SBM"))
  expect_true(inherits(mySBM, "SimpleSBM"))
  expect_true(inherits(mySBM, "SimpleSBM_fit"))

  ## Checking field access and format prior to estimation
  ## parameters
  expect_equal(mySBM$modelName, 'gaussian')
  expect_equal(unname(mySBM$nbNodes), nbNodes)
  expect_equal(mySBM$dimLabels, c(node="nodeName"))
  expect_equal(mySBM$nbDyads, nbNodes*(nbNodes - 1)/2)
  expect_true(all(is.na(diag(mySBM$networkData))))
  expect_true(isSymmetric(mySBM$networkData))
  expect_true(!mySBM$directed)
  expect_true(is.matrix(mySBM$connectParam$mean))

  ## covariates
  expect_true(all(dim(mySBM$covarEffect) == c(nbNodes, nbNodes)))
  expect_equal(mySBM$nbCovariates, 2)
  expect_equal(mySBM$covarList, covarList)
  expect_equal(mySBM$covarParam, c(0,0))

  ## S3 methods
  expect_equal(coef(mySBM, 'connectivity'), mySBM$connectParam)
  expect_equal(coef(mySBM, 'block')       , mySBM$blockProp)
  expect_equal(coef(mySBM, 'covariates')  , mySBM$covarParam)

  ## Estimation-----------------------------------------------------------------
  BM_out <- mySBM$optimize(estimOptions=list(verbosity = 0))
  mySBM$setModel(2)

  ## Expectation
  expect_equal(dim(mySBM$expectation), c(nbNodes, nbNodes))
  expect_gt(mySBM$connectParam$var, 0)

  ## blocks
  expect_equal(mySBM$nbBlocks, nbBlocks)
  expect_equal(dim(mySBM$probMemberships), c(nbNodes, nbBlocks))
  expect_equal(sort(unique(mySBM$memberships)), 1:nbBlocks)
  expect_equal(length(mySBM$memberships), nbNodes)

  ## correctness
  expect_lt(rmse(mySBM$connectParam$mean, means), 1e-1)
  expect_lt(rmse(mySBM$covarParam, covarParam), 0.1)
  expect_lt(1 - aricode::ARI(mySBM$memberships, mySampler$memberships), 1e-1)

  ## S3 methods
  expect_equal(coef(mySBM, 'connectivity'), mySBM$connectParam)
  expect_equal(coef(mySBM, 'block')       , mySBM$blockProp)
  expect_equal(coef(mySBM, 'covariates')  , mySBM$covarParam)
  expect_equal(mySBM$predict(), predict(mySBM))
  expect_equal(predict(mySBM), fitted(mySBM))

  expect_equal(predict(mySBM, covarList), fitted(mySBM))

  ## prediction wrt BM
  for (Q in mySBM$storedModels$indexModel) {
    pred_bm  <- BM_out$prediction(Q = Q)
    mySBM$setModel(Q)
    pred_sbm <- predict(mySBM)
    expect_lt( rmse(pred_bm, pred_sbm), 1e-12)
  }

})


test_that("SimpleSBM_fit 'Gaussian' model, undirected, two covariates", {

  ## SIMPLE UNDIRECTED GAUSSIAN SBM
  means <- matrix(c(1, 4, 7, 10),2,2)
  connectParam <- list(mean = means, var = 2)

  ## Basic construction - check for wrong specifications
  mySampler <- SimpleSBM$new('gaussian', nbNodes, TRUE, blockProp, connectParam, covarParam = covarParam, covarList = covarList_directed)
  mySampler$rMemberships(store = TRUE)
  mySampler$rEdges(store = TRUE)

  ## Construction----------------------------------------------------------------
  mySBM <- SimpleSBM_fit$new(mySampler$networkData, 'gaussian', TRUE, covarList = covarList_directed)
  expect_error(SimpleSBM_fit$new(SamplerBernoulli$networkData, 'normal', TRUE, covarList = covarList_directed))
  expect_error(SimpleSBM_fit$new(SamplerBernoulli$networkData[1:20, 1:30], 'gaussian', TRUE, covarList = covarList_directed))
  expect_error(SimpleSBM_fit$new(SamplerBernoulli$networkData, 'gaussian', FALSE, covarList = covarList_directed))

  ## Checking class
  expect_true(inherits(mySBM, "SBM"))
  expect_true(inherits(mySBM, "SimpleSBM"))
  expect_true(inherits(mySBM, "SimpleSBM_fit"))

  ## Checking field access and format prior to estimation
  ## parameters
  expect_equal(mySBM$modelName, 'gaussian')
  expect_equal(unname(mySBM$nbNodes), nbNodes)
  expect_equal(mySBM$dimLabels, c(node="nodeName"))
  expect_equal(mySBM$nbDyads, nbNodes*(nbNodes - 1))
  expect_true(all(is.na(diag(mySBM$networkData))))
  expect_true(!isSymmetric(mySBM$networkData))
  expect_true(mySBM$directed)
  expect_true(is.matrix(mySBM$connectParam$mean))

  ## covariates
  expect_true(all(dim(mySBM$covarEffect) == c(nbNodes, nbNodes)))
  expect_equal(mySBM$nbCovariates, 2)
  expect_equal(mySBM$covarList, covarList_directed)
  expect_equal(mySBM$covarParam, c(0,0))

  ## S3 methods
  expect_equal(coef(mySBM, 'connectivity'), mySBM$connectParam)
  expect_equal(coef(mySBM, 'block')       , mySBM$blockProp)
  expect_equal(coef(mySBM, 'covariates')  , mySBM$covarParam)

  ## Estimation-----------------------------------------------------------------
  BM_out <- mySBM$optimize(estimOptions=list(verbosity = 0))
  mySBM$setModel(2)

  ## Expectation
  expect_equal(dim(mySBM$expectation), c(nbNodes, nbNodes))
  expect_gt(mySBM$connectParam$var, 0)

  ## blocks
  expect_equal(mySBM$nbBlocks, nbBlocks)
  expect_equal(dim(mySBM$probMemberships), c(nbNodes, nbBlocks))
  expect_equal(sort(unique(mySBM$memberships)), 1:nbBlocks)
  expect_equal(length(mySBM$memberships), nbNodes)

  ## correctness
  expect_lt(rmse(sort(mySBM$connectParam$mean), means), 1e-1)
  expect_lt(rmse(mySBM$covarParam, covarParam), 0.1)
  expect_lt(1 - aricode::ARI(mySBM$memberships, mySampler$memberships), 1e-1)

  ## S3 methods
  expect_equal(coef(mySBM, 'connectivity'), mySBM$connectParam)
  expect_equal(coef(mySBM, 'block')       , mySBM$blockProp)
  expect_equal(coef(mySBM, 'covariates')  , mySBM$covarParam)
  expect_equal(mySBM$predict(), predict(mySBM))
  expect_equal(predict(mySBM), fitted(mySBM))

  expect_equal(predict(mySBM, covarList), fitted(mySBM))

  ## prediction wrt BM
  for (Q in mySBM$storedModels$indexModel) {
    pred_bm  <- BM_out$prediction(Q = Q)
    mySBM$setModel(Q)
    pred_sbm <- predict(mySBM)
    expect_lt( rmse(pred_bm, pred_sbm), 1e-12)
  }

})
