set.seed(1234)

rmse <- function(theta, theta_star) { sqrt(sum((theta - theta_star)^2)/sum(theta_star^2)) }

## Common parameters
nbNodes  <- 30
nbBlocks <- 2
blockProp <- c(.75, .25) # group proportions
covarParam <- c(-2,2)
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
  mySampler <- SimpleSBM_sampler$new('bernoulli', nbNodes, FALSE, blockProp, connectParam, covarParam[1], covarList[1])

  ## Construction----------------------------------------------------------------
  mySBM <- SimpleSBM_fit$new(mySampler$netMatrix, 'bernoulli', FALSE, covarList[1])
  expect_error(SimpleSBM_fit$new(mySampler$netMatrix, 'bernouilli', FALSE, covarList[1]))
  expect_error(SimpleSBM_fit$new(mySampler$netMatrix[1:20, 1:30], 'bernouilli', FALSE, covarList[1]))
  expect_error(SimpleSBM_fit$new(mySampler$netMatrix, 'bernoulli', TRUE, covarList[1]))
  expect_error(SimpleSBM_fit$new(mySampler$netMatrix, 'bernoulli', FALSE, covarList[[1]]))

  ## Checking class
  expect_true(inherits(mySBM, "SBM"))
  expect_true(inherits(mySBM, "SBM_fit"))
  expect_true(inherits(mySBM, "SimpleSBM_fit"))

  ## Checking field access and format prior to estimation
  ## parameters
  expect_equal(mySBM$modelName, 'bernoulli')
  expect_equal(mySBM$nbNodes, nbNodes)
  expect_equal(mySBM$dimension, c(nbNodes, nbNodes))
  expect_equal(mySBM$nbDyads, nbNodes*(nbNodes - 1)/2)
  expect_true(all(is.na(diag(mySBM$netMatrix))))
  expect_true(isSymmetric(mySBM$netMatrix))
  expect_true(!mySBM$directed)
  expect_true(is.na(mySBM$connectParam$mean))

  ## covariates
  expect_null(mySBM$covarExpect)
  expect_equal(mySBM$nbCovariates, 1)
  expect_equal(mySBM$covarList, covarList[1])
  expect_equal(mySBM$covarParam, c(0))

  ## Estimation-----------------------------------------------------------------
  mySBM$optimize(verbosity = 0, fast = TRUE)
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

})

test_that("SimpleSBM_fit 'Bernoulli' model, directed, one covariate", {

  ## SIMPLE UNDIRECTED BERNOULLI SBM
  means <- matrix(c(0.1, 0.4, 0.6, 0.9), 2,  2)
  connectParam <- list(mean = means)

  ## Basic construction - check for wrong specifications
  mySampler <- SimpleSBM_sampler$new('bernoulli', nbNodes, TRUE, blockProp, connectParam, covarParam[1], covarList_directed[1])

  ## Construction----------------------------------------------------------------
  mySBM <- SimpleSBM_fit$new(mySampler$netMatrix, 'bernoulli', TRUE, covarList_directed[1])
  expect_error(SimpleSBM_fit$new(mySampler$netMatrix, 'bernouilli', TRUE, covarList_directed[1]))
  expect_error(SimpleSBM_fit$new(mySampler$netMatrix, 'bernoulli', FALSE, covarList_directed[1]))
  expect_error(SimpleSBM_fit$new(mySampler$netMatrix[1:20, 1:30], 'bernouilli', FALSE, covarList_directed[1]))

  ## Checking class
  expect_true(inherits(mySBM, "SBM"))
  expect_true(inherits(mySBM, "SBM_fit"))
  expect_true(inherits(mySBM, "SimpleSBM_fit"))

  ## Checking field access and format prior to estimation
  ## parameters
  expect_equal(mySBM$modelName, 'bernoulli')
  expect_equal(mySBM$nbNodes, nbNodes)
  expect_equal(mySBM$dimension, c(nbNodes, nbNodes))
  expect_equal(mySBM$nbDyads, nbNodes*(nbNodes - 1))
  expect_true(all(is.na(diag(mySBM$netMatrix))))
  expect_true(!isSymmetric(mySBM$netMatrix))
  expect_true(mySBM$directed)
  expect_true(is.na(mySBM$connectParam$mean))

  ## covariates
  expect_null(mySBM$covarExpect)
  expect_equal(mySBM$nbCovariates, 1)
  expect_equal(mySBM$covarList, covarList_directed[1])
  expect_equal(mySBM$covarParam, c(0))

  ## Estimation-----------------------------------------------------------------
  mySBM$optimize(verbosity = 0, fast = TRUE)
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

})

test_that("SimpleSBM_fit 'Poisson' model, undirected, two covariates", {

  ## SIMPLE UNDIRECTED POISSON SBM
  means <- diag(15., 2) + 5
  connectParam <- list(mean = means)

  ## Basic construction - check for wrong specifications
  mySampler <- SimpleSBM_sampler$new('poisson', nbNodes, FALSE, blockProp, connectParam, covarParam, covarList)

  ## Construction----------------------------------------------------------------
  mySBM <- SimpleSBM_fit$new(mySampler$netMatrix, 'poisson', FALSE, covarList)
  expect_error(SimpleSBM_fit$new(SamplerBernoulli$netMatrix, 'poison', FALSE, covarList))
  expect_error(SimpleSBM_fit$new(SamplerBernoulli$netMatrix[1:20, 1:30], 'poisson', FALSE, covarList))
  expect_error(SimpleSBM_fit$new(SamplerBernoulli$netMatrix, 'poisson', TRUE, covarList))

  ## Checking class
  expect_true(inherits(mySBM, "SBM"))
  expect_true(inherits(mySBM, "SBM_fit"))
  expect_true(inherits(mySBM, "SimpleSBM_fit"))

  ## Checking field access and format prior to estimation
  ## parameters
  expect_equal(mySBM$modelName, 'poisson')
  expect_equal(mySBM$nbNodes, nbNodes)
  expect_equal(mySBM$dimension, c(nbNodes, nbNodes))
  expect_equal(mySBM$nbDyads, nbNodes*(nbNodes - 1)/2)
  expect_true(all(is.na(diag(mySBM$netMatrix))))
  expect_true(isSymmetric(mySBM$netMatrix))
  expect_true(!mySBM$directed)
  expect_true(is.na(mySBM$connectParam$mean))

  ## covariates
  expect_null(mySBM$covarExpect)
  expect_equal(mySBM$nbCovariates, 2)
  expect_equal(mySBM$covarList, covarList)
  expect_equal(mySBM$covarParam, c(0,0))

  ## Estimation-----------------------------------------------------------------
  mySBM$optimize(verbosity = 0)
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

})

test_that("SimpleSBM_fit 'Poisson' model, directed, two covariates", {

  ## SIMPLE UNDIRECTED BERNOULLI SBM
  means <- matrix(c(1, 4, 7, 9), 2,  2)
  connectParam <- list(mean = means)

  ## Basic construction - check for wrong specifications
  mySampler <- SimpleSBM_sampler$new('poisson', nbNodes, TRUE, blockProp, connectParam, covarParam, covarList_directed)

  ## Construction----------------------------------------------------------------
  mySBM <- SimpleSBM_fit$new(mySampler$netMatrix, 'poisson', TRUE, covarList_directed)
  expect_error(SimpleSBM_fit$new(SamplerBernoulli$netMatrix, 'poison', TRUE, covarList_directed))
  expect_error(SimpleSBM_fit$new(SamplerBernoulli$netMatrix[1:20, 1:30], 'poisson', FALSE, covarList_directed))
  expect_error(SimpleSBM_fit$new(SamplerBernoulli$netMatrix, 'poisson', FALSE, covarList_directed))

  ## Checking class
  expect_true(inherits(mySBM, "SBM"))
  expect_true(inherits(mySBM, "SBM_fit"))
  expect_true(inherits(mySBM, "SimpleSBM_fit"))

  ## Checking field access and format prior to estimation
  ## parameters
  expect_equal(mySBM$modelName, 'poisson')
  expect_equal(mySBM$nbNodes, nbNodes)
  expect_equal(mySBM$dimension, c(nbNodes, nbNodes))
  expect_equal(mySBM$nbDyads, nbNodes*(nbNodes - 1))
  expect_true(all(is.na(diag(mySBM$netMatrix))))
  expect_true(!isSymmetric(mySBM$netMatrix))
  expect_true(mySBM$directed)
  expect_true(is.na(mySBM$connectParam$mean))

  ## covariates
  expect_null(mySBM$covarExpect)
  expect_equal(mySBM$nbCovariates, 2)
  expect_equal(mySBM$covarList, covarList_directed)
  expect_equal(mySBM$covarParam, c(0,0))

  ## Estimation-----------------------------------------------------------------
  mySBM$optimize(verbosity = 0)
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

})


test_that("SimpleSBM_fit 'Gaussian' model, undirected, two covariates", {

  ## SIMPLE UNDIRECTED GAUSSIAN SBM
  means <- diag(15., 2) + 5 # connectivity matrix: affiliation network
  connectParam <- list(mean = means, var = 2)

  ## Basic construction - check for wrong specifications
  mySampler <- SimpleSBM_sampler$new('gaussian', nbNodes, FALSE, blockProp, connectParam, covarParam, covarList)

  ## Construction----------------------------------------------------------------
  mySBM <- SimpleSBM_fit$new(mySampler$netMatrix, 'gaussian', FALSE, covarList)
  expect_error(SimpleSBM_fit$new(SamplerBernoulli$netMatrix, 'normal', FALSE, covarList))
  expect_error(SimpleSBM_fit$new(SamplerBernoulli$netMatrix[1:20, 1:30], 'gaussian', FALSE, covarList))
  expect_error(SimpleSBM_fit$new(SamplerBernoulli$netMatrix, 'gaussian', TRUE, covarList))

  ## Checking class
  expect_true(inherits(mySBM, "SBM"))
  expect_true(inherits(mySBM, "SBM_fit"))
  expect_true(inherits(mySBM, "SimpleSBM_fit"))

  ## Checking field access and format prior to estimation
  ## parameters
  expect_equal(mySBM$modelName, 'gaussian')
  expect_equal(mySBM$nbNodes, nbNodes)
  expect_equal(mySBM$dimension, c(nbNodes, nbNodes))
  expect_equal(mySBM$nbDyads, nbNodes*(nbNodes - 1)/2)
  expect_true(all(is.na(diag(mySBM$netMatrix))))
  expect_true(isSymmetric(mySBM$netMatrix))
  expect_true(!mySBM$directed)
  expect_true(is.na(mySBM$connectParam$mean))

  ## covariates
  expect_null(mySBM$covarExpect)
  expect_equal(mySBM$nbCovariates, 2)
  expect_equal(mySBM$covarList, covarList)
  expect_equal(mySBM$covarParam, c(0,0))

  ## Estimation-----------------------------------------------------------------
  mySBM$optimize(verbosity = 0)
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

})





test_that("SimpleSBM_fit 'Gaussian' model, undirected, two covariates", {

  ## SIMPLE UNDIRECTED GAUSSIAN SBM
  means <- matrix(c(1, 4, 7, 10),2,2)
  connectParam <- list(mean = means, var = 2)

  ## Basic construction - check for wrong specifications
  mySampler <- SimpleSBM_sampler$new('gaussian', nbNodes, TRUE, blockProp, connectParam, covarParam, covarList_directed)

  ## Construction----------------------------------------------------------------
  mySBM <- SimpleSBM_fit$new(mySampler$netMatrix, 'gaussian', TRUE, covarList_directed)
  expect_error(SimpleSBM_fit$new(SamplerBernoulli$netMatrix, 'normal', TRUE, covarList_directed))
  expect_error(SimpleSBM_fit$new(SamplerBernoulli$netMatrix[1:20, 1:30], 'gaussian', TRUE, covarList_directed))
  expect_error(SimpleSBM_fit$new(SamplerBernoulli$netMatrix, 'gaussian', FALSE, covarList_directed))

  ## Checking class
  expect_true(inherits(mySBM, "SBM"))
  expect_true(inherits(mySBM, "SBM_fit"))
  expect_true(inherits(mySBM, "SimpleSBM_fit"))

  ## Checking field access and format prior to estimation
  ## parameters
  expect_equal(mySBM$modelName, 'gaussian')
  expect_equal(mySBM$nbNodes, nbNodes)
  expect_equal(mySBM$dimension, c(nbNodes, nbNodes))
  expect_equal(mySBM$nbDyads, nbNodes*(nbNodes - 1))
  expect_true(all(is.na(diag(mySBM$netMatrix))))
  expect_true(!isSymmetric(mySBM$netMatrix))
  expect_true(mySBM$directed)
  expect_true(is.na(mySBM$connectParam$mean))

  ## covariates
  expect_null(mySBM$covarExpect)
  expect_equal(mySBM$nbCovariates, 2)
  expect_equal(mySBM$covarList, covarList_directed)
  expect_equal(mySBM$covarParam, c(0,0))

  ## Estimation-----------------------------------------------------------------
  mySBM$optimize(verbosity = 0)
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

})