set.seed(1234)

## Common parameters
nbNodes  <- c(100, 120)
blockProp <- list(row = c(.5, .5), col = c(1/3, 1/3, 1/3)) # group proportions
nbBlocks <- sapply(blockProp, length)

test_that("Construction, fields access and other basics work in class BipartiteSBM_Sampler (Bernoulli model, no covariate)", {

  ## BIPARTITE UNDIRECTED BERNOULLI SBM
  means <- matrix(runif(6), 2, 3)  # connectivity matrix
  connectParam <- list(mean = means)

  ## Basic construction - check for wrong specifications
  mySampler <- BipartiteSBM$new('bernoulli', nbNodes, blockProp, connectParam)
  expect_error(SimpleSBM$new('bernouilli',nbNodes, blockProp, connectParam))
  expect_error(SimpleSBM$new('bernoulli', -10    , blockProp, connectParam))
  expect_error(SimpleSBM$new('bernoulli', c(1,2) , blockProp, connectParam))
  expect_error(SimpleSBM$new('bernoulli', nbNodes, 2, connectParam))
  expect_error(SimpleSBM$new('bernoulli', nbNodes, c(0,1), connectParam))
  expect_error(SimpleSBM$new('bernoulli', nbNodes, blockProp, list(mean = matrix( 2, nbBlocks[1], nbBlocks[2]))))
  expect_error(SimpleSBM$new('bernoulli', nbNodes, blockProp, list(mean = matrix(-2, nbBlocks[1], nbBlocks[2]))))
  expect_error(SimpleSBM$new('bernoulli', nbNodes, FALSE, blockProp, list(mean = matrix(0, nbBlocks[1] - 1, nbBlocks[2]))))

  ## Checking class
  expect_true(inherits(mySampler, "SBM"))
  expect_true(inherits(mySampler, "BipartiteSBM"))

  ## Checking field access and format

  ## parameters
  expect_equal(mySampler$modelName, 'bernoulli')
  expect_equal(unname(mySampler$nbNodes), nbNodes)
  expect_equal(mySampler$nbDyads, nbNodes[1]*nbNodes[2])

  ## covariates
  expect_null(mySampler$covarExpect)
  expect_equal(mySampler$nbCovariates, 0)
  expect_equal(mySampler$covarList, list())
  expect_equal(mySampler$covarParam, numeric(0))

  expect_equal(mySampler$connectParam$mean, means)
  expect_null(mySampler$connectParam$var)

  ## network
  mySampler$rMemberships(store = TRUE)
  mySampler$rEdges(store = TRUE)
  expect_equal(dim(mySampler$expectation), nbNodes)
  expect_true(all(mySampler$expectation >= 0, na.rm = TRUE))
  expect_true(all(mySampler$expectation <= 1, na.rm = TRUE))
  expect_false(isSymmetric(mySampler$networkData))

  ## blocks
  expect_equal(mySampler$blockProp, blockProp)
  expect_equal(mySampler$nbBlocks, nbBlocks)
  expect_equivalent(dim(mySampler$indMemberships[[1]]), c(nbNodes[1], nbBlocks[1]))
  expect_equivalent(dim(mySampler$indMemberships[[2]]), c(nbNodes[2], nbBlocks[2]))
  expect_equal(sort(unique(mySampler$memberships[[1]])), 1:nbBlocks[1])
  expect_equal(sort(unique(mySampler$memberships[[2]])), 1:nbBlocks[2])
  expect_equal(length(mySampler$memberships[[1]]), nbNodes[1])
  expect_equal(length(mySampler$memberships[[2]]), nbNodes[2])

  ## covariates
  expect_equal(mySampler$nbCovariates, 0)
  expect_equal(mySampler$covarList, list())
  expect_equal(mySampler$covarParam, numeric(0))

})


test_that("Construction, fields access and other basics work in class BipartiteSBM_Sampler (Poisson model, no covariate)", {

  ## SIMPLE UNDIRECTED POISSON SBM
  means <- matrix(rbinom(6, 30, 0.25), 2, 3)  # connectivity matrix
  connectParam <- list(mean = means)

  ## Basic construction - check for wrong specifications
  mySampler <- BipartiteSBM$new('poisson', nbNodes, blockProp, connectParam)
  expect_error(BipartiteSBM$new('poison' , nbNodes, blockProp, connectParam))
  expect_error(BipartiteSBM$new('poisson', -10    , blockProp, connectParam))
  expect_error(BipartiteSBM$new('poisson', nbNodes, -2, connectParam))
  expect_error(BipartiteSBM$new('poisson', nbNodes,  c(0,1), connectParam))
  expect_error(BipartiteSBM$new('poisson', nbNodes, blockProp, list(mean = matrix(-2, nbBlocks[1], nbBlocks[2]))))
  expect_error(BipartiteSBM$new('poisson', nbNodes, blockProp, list(mean = matrix(2 , nbBlocks[1] - 1, nbBlocks[2]))))

  ## Checking class
  expect_true(inherits(mySampler, "SBM"))
  expect_true(inherits(mySampler, "BipartiteSBM"))

  ## Checking field access and format

  ## parameters
  expect_equal(mySampler$modelName, 'poisson')
  expect_equal(unname(mySampler$nbNodes), nbNodes)
  expect_equal(mySampler$nbDyads, nbNodes[1]*nbNodes[2])
  expect_equal(mySampler$connectParam$mean, means)
  expect_null(mySampler$connectParam$var)

  ## network
  mySampler$rMemberships(store = TRUE)
  mySampler$rEdges(store = TRUE)
  expect_equal(dim(mySampler$expectation), nbNodes)
  expect_true(all(mySampler$expectation >= 0, na.rm = TRUE))
  expect_false(isSymmetric(mySampler$networkData))

  ## blocks
  expect_equal(mySampler$blockProp, blockProp)
  expect_equal(mySampler$nbBlocks, nbBlocks)
  expect_equivalent(dim(mySampler$indMemberships[[1]]), c(nbNodes[1], nbBlocks[1]))
  expect_equivalent(dim(mySampler$indMemberships[[2]]), c(nbNodes[2], nbBlocks[2]))
  expect_equal(sort(unique(mySampler$memberships[[1]])), 1:nbBlocks[1])
  expect_equal(sort(unique(mySampler$memberships[[2]])), 1:nbBlocks[2])
  expect_equal(length(mySampler$memberships[[1]]), nbNodes[1])
  expect_equal(length(mySampler$memberships[[2]]), nbNodes[2])

  ## covariates
  expect_equal(mySampler$nbCovariates, 0)
  expect_equal(mySampler$covarList, list())
  expect_equal(mySampler$covarParam, numeric(0))
})

test_that("Construction, fields access and other basics work in class BipartiteSBM_Sampler (Gaussian model, no covariate)", {

  ## SIMPLE UNDIRECTED GAUSSIAN SBM
  means <- matrix(c(0.05, 0.95, 0.4, 0.98, 0.15, 0.6), 2, 3)  # connectivity matrix
  connectParam <- list(mean = means, var = .1)

  ## Basic construction - check for wrong specifications
  mySampler <- BipartiteSBM$new('gaussian', nbNodes, blockProp, connectParam)
  expect_error(BipartiteSBM$new('normal'  , nbNodes, blockProp, connectParam))
  expect_error(BipartiteSBM$new('gaussian', -10    , blockProp, connectParam))
  expect_error(BipartiteSBM$new('gaussian', nbNodes, -2, connectParam))
  expect_error(BipartiteSBM$new('gaussian', nbNodes,  c(0,1), connectParam))
  expect_error(BipartiteSBM$new('gaussian', nbNodes, blockProp, list(var = -1, mean = means)))
  expect_error(BipartiteSBM$new('gaussian', nbNodes, blockProp, list(mean = matrix(runif(nbBlocks**2), nbBlocks, nbBlocks))))
  expect_error(BipartiteSBM$new('gaussian', nbNodes, blockProp, list(var = 1 , mean = matrix(2 , nbBlocks - 1, nbBlocks))))
  expect_error(BipartiteSBM$new('gaussian', nbNodes, blockProp, list(mean = matrix(-2, nbBlocks[1], nbBlocks[2]))))
  expect_error(BipartiteSBM$new('gaussian', nbNodes, blockProp, list(mean = matrix(2 , nbBlocks[1] - 1, nbBlocks[2]))))

  ## Checking class
  expect_true(inherits(mySampler, "SBM"))
  expect_true(inherits(mySampler, "BipartiteSBM"))

  ## Checking field access and format

  ## parameters
  expect_equal(mySampler$modelName, 'gaussian')
  expect_equal(unname(mySampler$nbNodes), nbNodes)
  expect_equal(mySampler$nbDyads, nbNodes[1]*nbNodes[2])
  expect_equal(mySampler$connectParam$mean, means)
  expect_gt(mySampler$connectParam$var, 0)

  ## network
  mySampler$rMemberships(store = TRUE)
  mySampler$rEdges(store = TRUE)
  expect_equal(dim(mySampler$expectation), nbNodes)
  expect_false(isSymmetric(mySampler$networkData))

  ## blocks
  expect_equal(mySampler$blockProp, blockProp)
  expect_equal(mySampler$nbBlocks, nbBlocks)
  expect_equivalent(dim(mySampler$indMemberships[[1]]), c(nbNodes[1], nbBlocks[1]))
  expect_equivalent(dim(mySampler$indMemberships[[2]]), c(nbNodes[2], nbBlocks[2]))
  expect_equal(sort(unique(mySampler$memberships[[1]])), 1:nbBlocks[1])
  expect_equal(sort(unique(mySampler$memberships[[2]])), 1:nbBlocks[2])
  expect_equal(length(mySampler$memberships[[1]]), nbNodes[1])
  expect_equal(length(mySampler$memberships[[2]]), nbNodes[2])

  ## covariates
  expect_equal(mySampler$nbCovariates, 0)
  expect_equal(mySampler$covarList, list())
  expect_equal(mySampler$covarParam, numeric(0))
})

