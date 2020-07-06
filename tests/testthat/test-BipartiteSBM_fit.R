rmse <- function(theta, theta_star) { sqrt(sum((theta - theta_star)^2)/sum(theta_star^2)) }

## Common parameters
nbNodes  <- c(100, 120)
blockProp <- list(row = c(.5, .5), col = c(1/3, 1/3, 1/3)) # group proportions
nbBlocks <- sapply(blockProp, length)

test_that("BipartiteSBM_fit 'Bernoulli' model, undirected, no covariate", {

  ## BIPARTITE UNDIRECTED BERNOULLI SBM
  means <- matrix(c(0.05, 0.95, 0.4, 0.75, 0.15, 0.6), 2, 3)  # connectivity matrix
  connectParam <- list(mean = means)

  ## Basic construction - check for wrong specifications
  mySampler <- BipartiteSBM_sampler$new('bernoulli', nbNodes, blockProp, connectParam)

  ## Construction----------------------------------------------------------------
  mySBM <- BipartiteSBM_fit$new(mySampler$netMatrix, 'bernoulli')
  expect_error(BipartiteSBM_fit$new(SamplerBernoulli$netMatrix, 'bernouilli'))

  ## Checking class
  expect_true(inherits(mySBM, "SBM"))
  expect_true(inherits(mySBM, "SBM_fit"))
  expect_true(inherits(mySBM, "BipartiteSBM_fit"))

  ## Checking field access and format prior to estimation
  ## parameters
  expect_equal(mySBM$modelName, 'bernoulli')
  expect_equal(mySBM$nbNodes, nbNodes)
  expect_equal(mySBM$dimension, nbNodes)
  expect_equal(mySBM$nbDyads, nbNodes[1]*nbNodes[2])
  expect_true(is.na(mySBM$connectParam$mean))

  ## covariates
  expect_null(mySBM$covarExpect)
  expect_equal(mySBM$nbCovariates, 0)
  expect_equal(mySBM$covarList, list())
  expect_equal(mySBM$covarParam, numeric(0))

  ## Estimation-----------------------------------------------------------------
  mySBM$optimize(verbosity = 0)
  mySBM$setModel(5)

  ## Expectation
  expect_equal(dim(mySBM$expectation), nbNodes)
  expect_true(all(mySBM$expectation >= 0, na.rm = TRUE))
  expect_true(all(mySBM$expectation <= 1, na.rm = TRUE))
  expect_null(mySBM$connectParam$var)

  ## blocks
  expect_equal(mySBM$nbBlocks, nbBlocks)
  expect_equivalent(dim(mySBM$probMemberships[[1]]), c(nbNodes[1], nbBlocks[1]))
  expect_equivalent(dim(mySBM$probMemberships[[2]]), c(nbNodes[2], nbBlocks[2]))
  expect_equal(sort(unique(mySBM$memberships[[1]])), 1:nbBlocks[1])
  expect_equal(sort(unique(mySBM$memberships[[2]])), 1:nbBlocks[2])

  ## correctness
  expect_lt(rmse(sort(mySBM$connectParam$mean), sort(means)), 1e-1)
  expect_lt(1 - aricode::ARI(mySBM$memberships[[1]], mySampler$memberships[[1]]), .1)
  expect_lt(1 - aricode::ARI(mySBM$memberships[[2]], mySampler$memberships[[2]]), .1)

})

test_that("BipartiteSBM_fit 'Poisson' model, undirected, no covariate", {

  ## SIMPLE DIRECTED POISSON SBM
  means <- matrix(c(10, 5, 7, 15, 20, 8), 2, 3)  # connectivity matrix
  connectParam <- list(mean = means)

  ## Basic construction - check for wrong specifications
  mySampler <- BipartiteSBM_sampler$new('poisson', nbNodes, blockProp, connectParam)

  ## Construction----------------------------------------------------------------
  mySBM <- BipartiteSBM_fit$new(mySampler$netMatrix, 'poisson')
  expect_error(BipartiteSBM_fit$new(SamplerBernoulli$netMatrix, 'poison'))

  ## Checking class
  expect_true(inherits(mySBM, "SBM"))
  expect_true(inherits(mySBM, "SBM_fit"))
  expect_true(inherits(mySBM, "BipartiteSBM_fit"))

  ## Checking field access and format prior to estimation
  ## parameters
  expect_equal(mySBM$modelName, 'poisson')
  expect_equal(mySBM$nbNodes, nbNodes)
  expect_equal(mySBM$dimension, nbNodes)
  expect_equal(mySBM$nbDyads, nbNodes[1]*nbNodes[2])
  expect_true(is.na(mySBM$connectParam$mean))

  ## covariates
  expect_null(mySBM$covarExpect)
  expect_equal(mySBM$nbCovariates, 0)
  expect_equal(mySBM$covarList, list())
  expect_equal(mySBM$covarParam, numeric(0))

  ## Estimation-----------------------------------------------------------------
  mySBM$optimize(verbosity = 0)
  mySBM$setModel(5)

  ## Expectation
  expect_equal(dim(mySBM$expectation), nbNodes)
  expect_true(all(mySBM$expectation >= 0, na.rm = TRUE))
  expect_null(mySBM$connectParam$var)

  ## blocks
  expect_equal(mySBM$nbBlocks, nbBlocks)
  expect_equivalent(dim(mySBM$probMemberships[[1]]), c(nbNodes[1], nbBlocks[1]))
  expect_equivalent(dim(mySBM$probMemberships[[2]]), c(nbNodes[2], nbBlocks[2]))
  expect_equal(sort(unique(mySBM$memberships[[1]])), 1:nbBlocks[1])
  expect_equal(sort(unique(mySBM$memberships[[2]])), 1:nbBlocks[2])

  ## correctness
  expect_lt(rmse(sort(mySBM$connectParam$mean), sort(means)), 1e-1)
  expect_lt(1 - aricode::ARI(mySBM$memberships[[1]], mySampler$memberships[[1]]), 1e-1)
  expect_lt(1 - aricode::ARI(mySBM$memberships[[2]], mySampler$memberships[[2]]), 1e-1)

})

test_that("BipartiteSBM_fit 'Gaussian' model, undirected, no covariate", {

  ## SIMPLE UNDIRECTED GAUSSIAN SBM
  means <- matrix(c(0.05, 0.95, 0.4, 0.98, 0.15, 0.6), 2, 3)  # connectivity matrix
  connectParam <- list(mean = means, var = .1)

  ## Basic construction - check for wrong specifications
  mySampler <- BipartiteSBM_sampler$new('gaussian', nbNodes, blockProp, connectParam)

  ## Construction----------------------------------------------------------------
  mySBM <- BipartiteSBM_fit$new(mySampler$netMatrix, 'gaussian')
  expect_error(BipartiteSBM_fit$new(SamplerBernoulli$netMatrix, 'groß'))

  ## Checking class
  expect_true(inherits(mySBM, "SBM"))
  expect_true(inherits(mySBM, "SBM_fit"))
  expect_true(inherits(mySBM, "BipartiteSBM_fit"))

  ## Checking field access and format prior to estimation
  ## parameters
  expect_equal(mySBM$modelName, 'gaussian')
  expect_equal(mySBM$nbNodes, nbNodes)
  expect_equal(mySBM$dimension, nbNodes)
  expect_equal(mySBM$nbDyads, nbNodes[1]*nbNodes[2])
  expect_true(is.na(mySBM$connectParam$mean))

  ## covariates
  expect_null(mySBM$covarExpect)
  expect_equal(mySBM$nbCovariates, 0)
  expect_equal(mySBM$covarList, list())
  expect_equal(mySBM$covarParam, numeric(0))

  ## Estimation-----------------------------------------------------------------
  mySBM$optimize(verbosity = 0)
  mySBM$setModel(5)

  ## Expectation
  expect_equal(dim(mySBM$expectation), nbNodes)
  expect_gt(mySBM$connectParam$var, 0)

  ## blocks
  expect_equal(mySBM$nbBlocks, nbBlocks)
  expect_equivalent(dim(mySBM$probMemberships[[1]]), c(nbNodes[1], nbBlocks[1]))
  expect_equivalent(dim(mySBM$probMemberships[[2]]), c(nbNodes[2], nbBlocks[2]))
  expect_equal(sort(unique(mySBM$memberships[[1]])), 1:nbBlocks[1])
  expect_equal(sort(unique(mySBM$memberships[[2]])), 1:nbBlocks[2])

  ## correctness
  expect_lt(rmse(sort(mySBM$connectParam$mean), sort(means)), 1e-1)
  expect_lt(1 - aricode::ARI(mySBM$memberships[[1]], mySampler$memberships[[1]]), 1e-1)
  expect_lt(1 - aricode::ARI(mySBM$memberships[[2]], mySampler$memberships[[2]]), 1e-1)

})

