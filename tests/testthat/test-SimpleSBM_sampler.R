## Common parameters
nbNodes  <- 90
nbBlocks <- 3
blockProp <- c(.5, .25, .25) # group proportions

test_that("Construction, fields access and other basics work in class SimpleSBM_Sampler (undirected Bernoulli model, no covariate)", {

  ## SIMPLE UNDIRECTED BERNOULLI SBM
  means <- diag(.4, 3) + 0.05
  connectParam <- list(mean = means)

  ## Basic construction - check for wrong specifications
  mySampler <- SimpleSBM_sampler$new('bernoulli', nbNodes, FALSE, blockProp, connectParam)
  expect_error(SimpleSBM_sampler$new('bernouilli',nbNodes, FALSE, blockProp, connectParam))
  expect_error(SimpleSBM_sampler$new('bernoulli', -10    , FALSE, blockProp, connectParam))
  expect_error(SimpleSBM_sampler$new('bernoulli', c(1,2) , FALSE, blockProp, connectParam))
  expect_error(SimpleSBM_sampler$new('bernoulli', nbNodes, FALSE, -2, connectParam))
  expect_error(SimpleSBM_sampler$new('bernoulli', nbNodes, FALSE,  c(0,1), connectParam))
  expect_error(SimpleSBM_sampler$new('bernoulli', nbNodes, FALSE, blockProp, list(mean = matrix( 2, nbBlocks, nbBlocks))))
  expect_error(SimpleSBM_sampler$new('bernoulli', nbNodes, FALSE, blockProp, list(mean = matrix(-2, nbBlocks, nbBlocks))))
  expect_error(SimpleSBM_sampler$new('bernoulli', nbNodes, FALSE, blockProp, list(mean = matrix(runif(nbBlocks**2), nbBlocks, nbBlocks))))
  expect_error(SimpleSBM_sampler$new('bernoulli', nbNodes, FALSE, blockProp, list(mean = matrix(0, nbBlocks - 1, nbBlocks))))

  ## Checking class
  expect_true(inherits(mySampler, "SBM"))
  expect_true(inherits(mySampler, "SBM_sampler"))
  expect_true(inherits(mySampler, "SimpleSBM_sampler"))

  ## Checking field access and format

  ## parameters
  expect_equal(mySampler$modelName, 'bernoulli')
  expect_equal(mySampler$nbNodes, nbNodes)
  expect_equal(mySampler$dimension, c(nbNodes, nbNodes))
  expect_equal(mySampler$nbDyads, nbNodes*(nbNodes - 1)/2)
  expect_equal(mySampler$connectParam$mean, means)
  expect_null(mySampler$connectParam$var)
  expect_equal(dim(mySampler$expectation), c(nbNodes, nbNodes))
  expect_true(all(mySampler$expectation >= 0, na.rm = TRUE))
  expect_true(all(mySampler$expectation <= 1, na.rm = TRUE))
  expect_true(isSymmetric(mySampler$netMatrix))
  expect_true(all(is.na(diag(mySampler$netMatrix))))
  expect_true(!mySampler$directed)

  ## blocks
  expect_equal(mySampler$blockProp, blockProp)
  expect_equal(mySampler$nbBlocks, nbBlocks)
  expect_equal(dim(mySampler$indMemberships), c(nbNodes, nbBlocks))
  expect_equal(sort(unique(mySampler$memberships)), 1:nbBlocks)
  expect_equal(length(mySampler$memberships), nbNodes)

  ## covariates
  expect_equal(mySampler$nbCovariates, 0)
  expect_equal(mySampler$covarList, list())
  expect_equal(mySampler$covarParam, numeric(0))

})

test_that("Construction, fields access and other basics work in class SimpleSBM_Sampler (directed Bernoulli model, no covariate)", {

  ## SIMPLE DIRECTED BERNOULLI SBM
  means <- matrix(runif(nbBlocks**2), nbBlocks, nbBlocks)
  connectParam <- list(mean = means)

  ## Basic construction - check for wrong specifications
  mySampler <- SimpleSBM_sampler$new('bernoulli', nbNodes, TRUE, blockProp, connectParam)
  expect_error(SimpleSBM_sampler$new('bernouilli',nbNodes, TRUE, blockProp, connectParam))
  expect_error(SimpleSBM_sampler$new('bernouilli',nbNodes, FALSE, blockProp, connectParam))
  expect_error(SimpleSBM_sampler$new('bernoulli', -10    , TRUE, blockProp, connectParam))
  expect_error(SimpleSBM_sampler$new('bernoulli', c(1,2) , TRUE, blockProp, connectParam))
  expect_error(SimpleSBM_sampler$new('bernoulli', nbNodes, TRUE, -2, connectParam))
  expect_error(SimpleSBM_sampler$new('bernoulli', nbNodes, TRUE,  c(0,1), connectParam))
  expect_error(SimpleSBM_sampler$new('bernoulli', nbNodes, TRUE, blockProp, list(mean = matrix( 2, nbBlocks, nbBlocks))))
  expect_error(SimpleSBM_sampler$new('bernoulli', nbNodes, TRUE, blockProp, list(mean = matrix(-2, nbBlocks, nbBlocks))))
  expect_error(SimpleSBM_sampler$new('bernoulli', nbNodes, TRUE, blockProp, list(mean = matrix(0, nbBlocks - 1, nbBlocks))))

  ## Checking class
  expect_true(inherits(mySampler, "SBM"))
  expect_true(inherits(mySampler, "SBM_sampler"))
  expect_true(inherits(mySampler, "SimpleSBM_sampler"))

  ## Checking field access and format

  ## parameters
  expect_equal(mySampler$modelName, 'bernoulli')
  expect_equal(mySampler$nbNodes, nbNodes)
  expect_equal(mySampler$dimension, c(nbNodes, nbNodes))
  expect_equal(mySampler$nbDyads, nbNodes*(nbNodes - 1))
  expect_equal(mySampler$connectParam$mean, means)
  expect_null(mySampler$connectParam$var)
  expect_equal(dim(mySampler$expectation), c(nbNodes, nbNodes))
  expect_true(all(mySampler$expectation >= 0, na.rm = TRUE))
  expect_true(all(mySampler$expectation <= 1, na.rm = TRUE))
  expect_true(all(is.na(diag(mySampler$netMatrix))))
  expect_true(!isSymmetric(mySampler$netMatrix))
  expect_true(mySampler$directed)

  ## blocks
  expect_equal(mySampler$blockProp, blockProp)
  expect_equal(mySampler$nbBlocks, nbBlocks)
  expect_equal(dim(mySampler$indMemberships), c(nbNodes, nbBlocks))
  expect_equal(sort(unique(mySampler$memberships)), 1:nbBlocks)
  expect_equal(length(mySampler$memberships), nbNodes)

  ## covariates
  expect_equal(mySampler$nbCovariates, 0)
  expect_equal(mySampler$covarList, list())
  expect_equal(mySampler$covarParam, numeric(0))

})

test_that("Construction, fields access and other basics work in class SimpleSBM_Sampler (undirected Poisson model, no covariate)", {

  ## SIMPLE UNDIRECTED POISSON SBM
  means <- diag(15., 3) + 5
  connectParam <- list(mean = means)

  ## Basic construction - check for wrong specifications
  mySampler <- SimpleSBM_sampler$new('poisson', nbNodes, FALSE, blockProp, connectParam)
  expect_error(SimpleSBM_sampler$new('poison' , nbNodes, FALSE, blockProp, connectParam))
  expect_error(SimpleSBM_sampler$new('poisson', -10    , FALSE, blockProp, connectParam))
  expect_error(SimpleSBM_sampler$new('poisson', c(1,2) , FALSE, blockProp, connectParam))
  expect_error(SimpleSBM_sampler$new('poisson', nbNodes, FALSE, -2, connectParam))
  expect_error(SimpleSBM_sampler$new('poisson', nbNodes, FALSE,  c(0,1), connectParam))
  expect_error(SimpleSBM_sampler$new('poisson', nbNodes, FALSE, blockProp, list(mean = matrix(-2, nbBlocks, nbBlocks))))
  expect_error(SimpleSBM_sampler$new('poisson', nbNodes, FALSE, blockProp, list(mean = round(40 * matrix(runif(nbBlocks**2), nbBlocks, nbBlocks)))))
  expect_error(SimpleSBM_sampler$new('poisson', nbNodes, FALSE, blockProp, list(mean = matrix(2 , nbBlocks - 1, nbBlocks))))

  ## Checking class
  expect_true(inherits(mySampler, "SBM"))
  expect_true(inherits(mySampler, "SBM_sampler"))
  expect_true(inherits(mySampler, "SimpleSBM_sampler"))

  ## Checking field access and format

  ## parameters
  expect_equal(mySampler$modelName, 'poisson')
  expect_equal(mySampler$nbNodes, nbNodes)
  expect_equal(mySampler$dimension, c(nbNodes, nbNodes))
  expect_equal(mySampler$nbDyads, nbNodes*(nbNodes - 1)/2)
  expect_equal(mySampler$connectParam$mean, means)
  expect_null(mySampler$connectParam$var)
  expect_equal(dim(mySampler$expectation), c(nbNodes, nbNodes))
  expect_true(all(mySampler$expectation >= 0, na.rm = TRUE))
  expect_true(all(is.na(diag(mySampler$netMatrix))))
  expect_true(isSymmetric(mySampler$netMatrix))

  ## blocks
  expect_equal(mySampler$blockProp, blockProp)
  expect_equal(mySampler$nbBlocks, nbBlocks)
  expect_equal(dim(mySampler$indMemberships), c(nbNodes, nbBlocks))
  expect_equal(sort(unique(mySampler$memberships)), 1:nbBlocks)
  expect_equal(length(mySampler$memberships), nbNodes)

  ## covariates
  expect_equal(mySampler$nbCovariates, 0)
  expect_equal(mySampler$covarList, list())
  expect_equal(mySampler$covarParam, numeric(0))

})


test_that("Construction, fields access and other basics work in class SimpleSBM_Sampler (directed Poisson model, no covariate)", {

  ## SIMPLE DIRECTED POISSON SBM
  means <- round(40 * matrix(runif(nbBlocks**2), nbBlocks, nbBlocks))
  connectParam <- list(mean = means)

  ## Basic construction - check for wrong specifications
  mySampler <- SimpleSBM_sampler$new('poisson', nbNodes, TRUE, blockProp, connectParam)
  expect_error(SimpleSBM_sampler$new('poison' , nbNodes, TRUE, blockProp, connectParam))
  expect_error(SimpleSBM_sampler$new('poisson', -10    , TRUE, blockProp, connectParam))
  expect_error(SimpleSBM_sampler$new('poisson', c(1,2) , TRUE, blockProp, connectParam))
  expect_error(SimpleSBM_sampler$new('poisson', nbNodes, TRUE, -2, connectParam))
  expect_error(SimpleSBM_sampler$new('poisson', nbNodes, TRUE,  c(0,1), connectParam))
  expect_error(SimpleSBM_sampler$new('poisson', nbNodes, TRUE, blockProp, list(mean = matrix(-2, nbBlocks, nbBlocks))))
  expect_error(SimpleSBM_sampler$new('poisson', nbNodes, TRUE, blockProp, list(mean = matrix(2 , nbBlocks - 1, nbBlocks))))

  ## Checking class
  expect_true(inherits(mySampler, "SBM"))
  expect_true(inherits(mySampler, "SBM_sampler"))
  expect_true(inherits(mySampler, "SimpleSBM_sampler"))

  ## Checking field access and format

  ## parameters
  expect_equal(mySampler$modelName, 'poisson')
  expect_equal(mySampler$nbNodes, nbNodes)
  expect_equal(mySampler$dimension, c(nbNodes, nbNodes))
  expect_equal(mySampler$nbDyads, nbNodes*(nbNodes - 1))
  expect_equal(mySampler$connectParam$mean, means)
  expect_null(mySampler$connectParam$var)
  expect_equal(dim(mySampler$expectation), c(nbNodes, nbNodes))
  expect_true(all(mySampler$expectation >= 0, na.rm = TRUE))
  expect_true(all(is.na(diag(mySampler$netMatrix))))
  expect_true(!isSymmetric(mySampler$netMatrix))
  expect_true(mySampler$directed)

  ## blocks
  expect_equal(mySampler$blockProp, blockProp)
  expect_equal(mySampler$nbBlocks, nbBlocks)
  expect_equal(dim(mySampler$indMemberships), c(nbNodes, nbBlocks))
  expect_equal(sort(unique(mySampler$memberships)), 1:nbBlocks)
  expect_equal(length(mySampler$memberships), nbNodes)

  ## covariates
  expect_equal(mySampler$nbCovariates, 0)
  expect_equal(mySampler$covarList, list())
  expect_equal(mySampler$covarParam, numeric(0))

})


test_that("Construction, fields access and other basics work in class SimpleSBM_Sampler (undirected Gaussian model, no covariate)", {

  ## SIMPLE UNDIRECTED GAUSSIAN SBM
  means <- diag(15., 3) + 5 # connectivity matrix: affiliation network
  connectParam <- list(mean = means, var = 2)

  ## Basic construction - check for wrong specifications
  mySampler <- SimpleSBM_sampler$new('gaussian', nbNodes, FALSE, blockProp, connectParam)
  expect_error(SimpleSBM_sampler$new('normal'  , nbNodes, FALSE, blockProp, connectParam))
  expect_error(SimpleSBM_sampler$new('gaussian', -10    , FALSE, blockProp, connectParam))
  expect_error(SimpleSBM_sampler$new('gaussian', c(1,2) , FALSE, blockProp, connectParam))
  expect_error(SimpleSBM_sampler$new('gaussian', nbNodes, FALSE, -2, connectParam))
  expect_error(SimpleSBM_sampler$new('gaussian', nbNodes, FALSE,  c(0,1), connectParam))
  expect_error(SimpleSBM_sampler$new('gaussian', nbNodes, FALSE, blockProp, list(mean = means, var = -1)))
  expect_error(SimpleSBM_sampler$new('gaussian', nbNodes, FALSE, blockProp, list(mean = matrix(runif(nbBlocks**2), nbBlocks, nbBlocks))))
  expect_error(SimpleSBM_sampler$new('gaussian', nbNodes, FALSE, blockProp, list(mean = matrix(2 , nbBlocks - 1, nbBlocks), var = 1)))

  ## Checking class
  expect_true(inherits(mySampler, "SBM"))
  expect_true(inherits(mySampler, "SBM_sampler"))
  expect_true(inherits(mySampler, "SimpleSBM_sampler"))

  ## Checking field access and format

  ## parameters
  expect_equal(mySampler$modelName, 'gaussian')
  expect_equal(mySampler$nbNodes, nbNodes)
  expect_equal(mySampler$dimension, c(nbNodes, nbNodes))
  expect_equal(mySampler$nbDyads, nbNodes*(nbNodes - 1)/2)
  expect_equal(mySampler$connectParam$mean, means)
  expect_equal(mySampler$connectParam$var, 2)
  expect_equal(dim(mySampler$expectation), c(nbNodes, nbNodes))
  expect_true(all(is.na(diag(mySampler$netMatrix))))
  expect_true(isSymmetric(mySampler$netMatrix))

  ## blocks
  expect_equal(mySampler$blockProp, blockProp)
  expect_equal(mySampler$nbBlocks, nbBlocks)
  expect_equal(dim(mySampler$indMemberships), c(nbNodes, nbBlocks))
  expect_equal(sort(unique(mySampler$memberships)), 1:nbBlocks)
  expect_equal(length(mySampler$memberships), nbNodes)

  ## covariates
  expect_equal(mySampler$nbCovariates, 0)
  expect_equal(mySampler$covarList, list())
  expect_equal(mySampler$covarParam, numeric(0))

})

test_that("Construction, fields access and other basics work in class SimpleSBM_Sampler (directed Gaussian model, no covariate)", {

  ## SIMPLE DIRECTED GAUSSIAN SBM
  means <- matrix(runif(nbBlocks**2), nbBlocks, nbBlocks)
  connectParam <- list(mean = means, var = 2)

  ## Basic construction - check for wrong specifications
  mySampler <- SimpleSBM_sampler$new('gaussian', nbNodes, TRUE, blockProp, connectParam)
  expect_error(SimpleSBM_sampler$new('normal'  , nbNodes, TRUE, blockProp, connectParam))
  expect_error(SimpleSBM_sampler$new('gaussian', nbNodes, FALSE, blockProp, connectParam))
  expect_error(SimpleSBM_sampler$new('gaussian', -10    , TRUE, blockProp, connectParam))
  expect_error(SimpleSBM_sampler$new('gaussian', c(1,2) , TRUE, blockProp, connectParam))
  expect_error(SimpleSBM_sampler$new('gaussian', nbNodes, TRUE, -2, connectParam))
  expect_error(SimpleSBM_sampler$new('gaussian', nbNodes, TRUE,  c(0,1), connectParam))
  expect_error(SimpleSBM_sampler$new('gaussian', nbNodes, TRUE, blockProp, list(var = -1, mean = means)))
  expect_error(SimpleSBM_sampler$new('gaussian', nbNodes, TRUE, blockProp, list(var = 1, mean = matrix(2 , nbBlocks - 1, nbBlocks))))

  ## Checking class
  expect_true(inherits(mySampler, "SBM"))
  expect_true(inherits(mySampler, "SBM_sampler"))
  expect_true(inherits(mySampler, "SimpleSBM_sampler"))

  ## Checking field access and format

  ## parameters
  expect_equal(mySampler$modelName, 'gaussian')
  expect_equal(mySampler$nbNodes, nbNodes)
  expect_equal(mySampler$dimension, c(nbNodes, nbNodes))
  expect_equal(mySampler$nbDyads, nbNodes*(nbNodes - 1))
  expect_equal(mySampler$connectParam$mean, means)
  expect_equal(mySampler$connectParam$var, 2)
  expect_equal(dim(mySampler$expectation), c(nbNodes, nbNodes))
  expect_true(all(is.na(diag(mySampler$netMatrix))))
  expect_true(!isSymmetric(mySampler$netMatrix))
  expect_true(mySampler$directed)

  ## blocks
  expect_equal(mySampler$blockProp, blockProp)
  expect_equal(mySampler$nbBlocks, nbBlocks)
  expect_equal(dim(mySampler$indMemberships), c(nbNodes, nbBlocks))
  expect_equal(sort(unique(mySampler$memberships)), 1:nbBlocks)
  expect_equal(length(mySampler$memberships), nbNodes)

  ## covariates
  expect_equal(mySampler$nbCovariates, 0)
  expect_equal(mySampler$covarList, list())
  expect_equal(mySampler$covarParam, numeric(0))

})

