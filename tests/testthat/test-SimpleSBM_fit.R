rmse <- function(theta, theta_star) { sqrt(sum((theta - theta_star)^2)/sum(theta_star^2)) }

## Common parameters
nbNodes  <- 90
nbBlocks <- 3
blockProp <- c(.5, .25, .25) # group proportions

test_that("SimpleSBM_fit 'Bernoulli' model, undirected, no covariate", {

  ## SIMPLE UNDIRECTED BERNOULLI SBM
  means <- diag(.4, 3) + 0.05
  connectParam <- list(mean = means)

  ## Basic construction - check for wrong specifications
  mySampler <- SimpleSBM_sampler$new('bernoulli', nbNodes, FALSE, blockProp, connectParam)

  ## Construction----------------------------------------------------------------
  mySBM <- SimpleSBM_fit$new(mySampler$netMatrix, 'bernoulli', FALSE)
  expect_error(SimpleSBM_fit$new(SamplerBernoulli$netMatrix, 'bernouilli', FALSE))
  expect_error(SimpleSBM_fit$new(SamplerBernoulli$netMatrix[1:20, 1:30], 'bernouilli', FALSE))
  expect_error(SimpleSBM_fit$new(SamplerBernoulli$netMatrix, 'bernoulli', TRUE))

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
  expect_equal(mySBM$nbCovariates, 0)
  expect_equal(mySBM$covarList, list())
  expect_equal(mySBM$covarParam, numeric(0))

  ## Estimation-----------------------------------------------------------------
  mySBM$optimize(verbosity = 0)
  mySBM$setModel(3)

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

  ## correctness
  expect_lt(rmse(mySBM$connectParam$mean, means), 0.2)
  expect_lt(1 - aricode::ARI(mySBM$memberships, mySampler$memberships), 0.2)

})

test_that("SimpleSBM_fit 'Poisson' model, undirected, no covariate", {

  ## SIMPLE UNDIRECTED POISSON SBM
  means <- diag(15., 3) + 5
  connectParam <- list(mean = means)

  ## Basic construction - check for wrong specifications
  mySampler <- SimpleSBM_sampler$new('poisson', nbNodes, FALSE, blockProp, connectParam)

  ## Construction----------------------------------------------------------------
  mySBM <- SimpleSBM_fit$new(mySampler$netMatrix, 'poisson', FALSE)
  expect_error(SimpleSBM_fit$new(SamplerBernoulli$netMatrix, 'poison', FALSE))
  expect_error(SimpleSBM_fit$new(SamplerBernoulli$netMatrix[1:20, 1:30], 'poisson', FALSE))
  expect_error(SimpleSBM_fit$new(SamplerBernoulli$netMatrix, 'poisson', TRUE))

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
  expect_equal(mySBM$nbCovariates, 0)
  expect_equal(mySBM$covarList, list())
  expect_equal(mySBM$covarParam, numeric(0))

  ## Estimation-----------------------------------------------------------------
  mySBM$optimize(verbosity = 0)
  mySBM$setModel(3)

  ## Expectation
  expect_equal(dim(mySBM$expectation), c(nbNodes, nbNodes))
  expect_true(all(mySampler$expectation >= 0, na.rm = TRUE))
  expect_null(mySBM$connectParam$var)

  ## blocks
  expect_equal(mySBM$nbBlocks, nbBlocks)
  expect_equal(dim(mySBM$probMemberships), c(nbNodes, nbBlocks))
  expect_equal(sort(unique(mySBM$memberships)), 1:nbBlocks)
  expect_equal(length(mySBM$memberships), nbNodes)

  ## correctness
  expect_lt(rmse(mySBM$connectParam$mean, means), 1e-1)
  expect_lt(1 - aricode::ARI(mySBM$memberships, mySampler$memberships), 1e-1)

})

test_that("SimpleSBM_fit 'Gaussian' model, undirected, no covariate", {

  ## SIMPLE UNDIRECTED GAUSSIAN SBM
  means <- diag(15., 3) + 5 # connectivity matrix: affiliation network
  connectParam <- list(mean = means, var = 2)

  ## Basic construction - check for wrong specifications
  mySampler <- SimpleSBM_sampler$new('gaussian', nbNodes, FALSE, blockProp, connectParam)

  ## Construction----------------------------------------------------------------
  mySBM <- SimpleSBM_fit$new(mySampler$netMatrix, 'gaussian', FALSE)
  expect_error(SimpleSBM_fit$new(SamplerBernoulli$netMatrix, 'normal', FALSE))
  expect_error(SimpleSBM_fit$new(SamplerBernoulli$netMatrix[1:20, 1:30], 'gaussian', FALSE))
  expect_error(SimpleSBM_fit$new(SamplerBernoulli$netMatrix, 'gaussian', TRUE))

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
  expect_equal(mySBM$nbCovariates, 0)
  expect_equal(mySBM$covarList, list())
  expect_equal(mySBM$covarParam, numeric(0))

  ## Estimation-----------------------------------------------------------------
  mySBM$optimize(verbosity = 0)
  mySBM$setModel(3)

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
  expect_lt(1 - aricode::ARI(mySBM$memberships, mySampler$memberships), 1e-1)

})

