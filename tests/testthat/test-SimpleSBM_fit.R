## Common parameters
nbNodes  <- 30
nbBlocks <- 2
blockProp <- c(.6, .4) # group proportions

## SIMPLE UNDIRECTED BERNOULLI SBM
means <- diag(.4, 2) + 0.05
connectParam <- list(mean = means)

## Basic construction - check for wrong specifications
SamplerBernoulli <- SimpleSBM_sampler$new('bernoulli', nbNodes, FALSE, blockProp, connectParam)

rmse <- function(theta, theta_star) {
  sqrt(sum((theta - theta_star)^2))
}

test_that("SimpleSBM_fit 'Bernoulli' model, without covariates", {

  ## Construction----------------------------------------------------------------
  mySBM <- SimpleSBM_fit$new(SamplerBernoulli$netMatrix, 'bernoulli', FALSE)
  expect_error(SimpleSBM_fit$new(SamplerBernoulli$netMatrix, 'bernouilli', FALSE))
  expect_error(SimpleSBM_fit$new(SamplerBernoulli$netMatrix[1:20, 1:30], 'bernouilli', FALSE))
  expect_error(SimpleSBM_fit$new(SamplerBernoulli$netMatrix, 'bernoulli', TRUE))

  ## Checking class
  expect_true(inherits(mySBM, "SBM"))
  expect_true(inherits(mySBM, "SBM_fit"))
  expect_true(inherits(mySBM, "SimpleSBM_fit"))

  ## Checking field access and format
  ## parameters
  expect_equal(mySBM$modelName, 'bernoulli')
  expect_equal(mySBM$nbNodes, nbNodes)
  expect_equal(mySBM$dimension, c(nbNodes, nbNodes))
  expect_equal(mySBM$nbDyads, nbNodes*(nbNodes - 1)/2)
  expect_equal(dim(mySBM$expectation), c(nbNodes, nbNodes))
  expect_true(all(is.na(diag(mySBM$netMatrix))))
  expect_true(isSymmetric(mySBM$netMatrix))
  expect_true(!mySBM$directed)

  ## Estimation-----------------------------------------------------------------
  mySBM$optimize(verbosity = 0)

  ## Checking field access and format
  ## parameters

  expect_lt(rmse(mySBM$connectParam$mean, means), 1e-1)
  expect_null(mySBM$connectParam$var)

  # expect_true(all(mySBM$expectation >= 0, na.rm = TRUE))
  # expect_true(all(mySBM$expectation <= 1, na.rm = TRUE))


})
