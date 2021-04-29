set.seed(1234)

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
  mySampler <- BipartiteSBM$new('bernoulli', nbNodes, blockProp, connectParam)
  mySampler$rMemberships(store = TRUE)
  mySampler$rEdges(store = TRUE)

  ## Construction----------------------------------------------------------------
  mySBM <- BipartiteSBM_fit$new(mySampler$networkData, 'bernoulli')
  expect_error(BipartiteSBM_fit$new(SamplerBernoulli$networkData, 'bernouilli'))

  ## Checking class
  expect_true(inherits(mySBM, "SBM"))
  expect_true(inherits(mySBM, "BipartiteSBM"))
  expect_true(inherits(mySBM, "BipartiteSBM_fit"))

  ## Checking field access and format prior to estimation
  ## parameters
  expect_equal(mySBM$modelName, 'bernoulli')
  expect_equal(unname(mySBM$nbNodes), nbNodes)
  expect_equal(mySBM$nbDyads, nbNodes[1]*nbNodes[2])
  expect_true(is.matrix(mySBM$connectParam$mean))

  ## covariates
  expect_equal(mySBM$covarEffect, numeric(0))
  expect_equal(mySBM$nbCovariates, 0)
  expect_equal(mySBM$covarList, list())
  expect_equal(mySBM$covarParam, numeric(0))

  ## S3 methods
  expect_equal(coef(mySBM, 'connectivity'), mySBM$connectParam)
  expect_equal(coef(mySBM, 'block')       , mySBM$blockProp)
  expect_equal(coef(mySBM, 'covariates')  , mySBM$covarParam)

  ## Estimation-----------------------------------------------------------------
  BM_out <- mySBM$optimize(estimOptions=list(verbosity = 0))
  mySBM$setModel(4)

  expect_equal(mySBM$nbConnectParam, unname(nbBlocks[1] * nbBlocks[2]))
  expect_equal(mySBM$penalty, nbBlocks[1] * nbBlocks[2] * log(nbNodes[1] * nbNodes[2]) +  (nbBlocks[1] - 1) * log(nbNodes[1]) + (nbBlocks[2] - 1) * log(nbNodes[2]))
  expect_equal(mySBM$entropy, -sum(mySBM$probMemberships[[1]] * log(mySBM$probMemberships[[1]]))
                              -sum(mySBM$probMemberships[[2]] * log(mySBM$probMemberships[[2]])))

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

  ## S3 methods
  expect_equal(coef(mySBM, 'connectivity'), mySBM$connectParam)
  expect_equal(coef(mySBM, 'block')       , mySBM$blockProp)
  expect_equal(coef(mySBM, 'covariates')  , mySBM$covarParam)
  expect_equal(mySBM$predict(), predict(mySBM))
  expect_equal(fitted(mySBM), predict(mySBM))

  ## correctness
  expect_lt(rmse(sort(mySBM$connectParam$mean), sort(means)), .2)
  expect_lt(1 - aricode::ARI(mySBM$memberships[[1]], mySampler$memberships[[1]]), .2)
  expect_lt(1 - aricode::ARI(mySBM$memberships[[2]], mySampler$memberships[[2]]), .2)

  ## prediction wrt BM
  for (Q in mySBM$storedModels$indexModel) {
    pred_bm  <- BM_out$prediction(Q = Q)
    mySBM$setModel(Q-1)
    pred_sbm <- predict(mySBM)
    expect_lt( rmse(pred_bm, pred_sbm), 1e-12)
  }

})

test_that("BipartiteSBM_fit 'Poisson' model, undirected, no covariate", {

  ## SIMPLE DIRECTED POISSON SBM
  means <- matrix(c(10, 5, 7, 15, 20, 8), 2, 3)  # connectivity matrix
  connectParam <- list(mean = means)

  ## Basic construction - check for wrong specifications
  mySampler <- BipartiteSBM$new('poisson', nbNodes, blockProp, connectParam)
  mySampler$rMemberships(store = TRUE)
  mySampler$rEdges(store = TRUE)

  ## Construction----------------------------------------------------------------
  mySBM <- BipartiteSBM_fit$new(mySampler$networkData, 'poisson')
  expect_error(BipartiteSBM_fit$new(SamplerBernoulli$networkData, 'poison'))

  ## Checking class
  expect_true(inherits(mySBM, "SBM"))
  expect_true(inherits(mySBM, "BipartiteSBM"))
  expect_true(inherits(mySBM, "BipartiteSBM_fit"))

  ## Checking field access and format prior to estimation
  ## parameters
  expect_equal(mySBM$modelName, 'poisson')
  expect_equal(unname(mySBM$nbNodes), nbNodes)
  expect_equal(mySBM$nbDyads, nbNodes[1]*nbNodes[2])
  expect_true(is.matrix(mySBM$connectParam$mean))

  ## covariates
  expect_equal(mySBM$covarEffect, numeric(0))
  expect_equal(mySBM$nbCovariates, 0)
  expect_equal(mySBM$covarList, list())
  expect_equal(mySBM$covarParam, numeric(0))

  ## S3 methods
  expect_equal(coef(mySBM, 'connectivity'), mySBM$connectParam)
  expect_equal(coef(mySBM, 'block')       , mySBM$blockProp)
  expect_equal(coef(mySBM, 'covariates')  , mySBM$covarParam)

  ## Estimation-----------------------------------------------------------------
  BM_out <- mySBM$optimize(estimOptions=list(verbosity = 0))
  mySBM$setModel(4)

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

  ## S3 methods
  expect_equal(coef(mySBM, 'connectivity'), mySBM$connectParam)
  expect_equal(coef(mySBM, 'block')       , mySBM$blockProp)
  expect_equal(coef(mySBM, 'covariates')  , mySBM$covarParam)
  expect_equal(mySBM$predict(), predict(mySBM))
  expect_equal(fitted(mySBM), predict(mySBM))

  ## correctness
  expect_lt(rmse(sort(mySBM$connectParam$mean), sort(means)), 1e-1)
  expect_lt(1 - aricode::ARI(mySBM$memberships[[1]], mySampler$memberships[[1]]), 1e-1)
  expect_lt(1 - aricode::ARI(mySBM$memberships[[2]], mySampler$memberships[[2]]), 1e-1)

  ## prediction wrt BM
  for (Q in mySBM$storedModels$indexModel) {
    pred_bm  <- BM_out$prediction(Q = Q)
    mySBM$setModel(Q-1)
    pred_sbm <- predict(mySBM)
    expect_lt( rmse(pred_bm, pred_sbm), 1e-12)
  }

})

test_that("BipartiteSBM_fit 'Gaussian' model, undirected, no covariate", {

  ## SIMPLE UNDIRECTED GAUSSIAN SBM
  means <- matrix(c(0.05, 0.95, 0.4, 0.98, 0.15, 0.6), 2, 3)  # connectivity matrix
  connectParam <- list(mean = means, var = .1)

  ## Basic construction - check for wrong specifications
  mySampler <- BipartiteSBM$new('gaussian', nbNodes, blockProp, connectParam)
  mySampler$rMemberships(store = TRUE)
  mySampler$rEdges(store = TRUE)

  ## Construction----------------------------------------------------------------
  mySBM <- BipartiteSBM_fit$new(mySampler$networkData, 'gaussian')
  expect_error(BipartiteSBM_fit$new(SamplerBernoulli$networkData, 'groß'))

  ## Checking class
  expect_true(inherits(mySBM, "SBM"))
  expect_true(inherits(mySBM, "BipartiteSBM"))
  expect_true(inherits(mySBM, "BipartiteSBM_fit"))

  ## Checking field access and format prior to estimation
  ## parameters
  expect_equal(mySBM$modelName, 'gaussian')
  expect_equal(unname(mySBM$nbNodes), nbNodes)
  expect_equal(mySBM$nbDyads, nbNodes[1]*nbNodes[2])
  expect_true(is.matrix(mySBM$connectParam$mean))

  ## covariates
  expect_equal(mySBM$covarEffect, numeric(0))
  expect_equal(mySBM$nbCovariates, 0)
  expect_equal(mySBM$covarList, list())
  expect_equal(mySBM$covarParam, numeric(0))

  ## S3 methods
  expect_equal(coef(mySBM, 'connectivity'), mySBM$connectParam)
  expect_equal(coef(mySBM, 'block')       , mySBM$blockProp)
  expect_equal(coef(mySBM, 'covariates')  , mySBM$covarParam)

  ## Estimation-----------------------------------------------------------------
  BM_out <- mySBM$optimize(estimOptions=list(verbosity = 0))
  mySBM$setModel(4)

  ## Expectation
  expect_equal(dim(mySBM$expectation), nbNodes)
  expect_gt(mySBM$connectParam$var, 0)

  ## blocks
  expect_equal(mySBM$nbBlocks, nbBlocks)
  expect_equivalent(dim(mySBM$probMemberships[[1]]), c(nbNodes[1], nbBlocks[1]))
  expect_equivalent(dim(mySBM$probMemberships[[2]]), c(nbNodes[2], nbBlocks[2]))
  expect_equal(sort(unique(mySBM$memberships[[1]])), 1:nbBlocks[1])
  expect_equal(sort(unique(mySBM$memberships[[2]])), 1:nbBlocks[2])

  ## S3 methods
  expect_equal(coef(mySBM, 'connectivity'), mySBM$connectParam)
  expect_equal(coef(mySBM, 'block')       , mySBM$blockProp)
  expect_equal(coef(mySBM, 'covariates')  , mySBM$covarParam)
  expect_equal(mySBM$predict(), predict(mySBM))
  expect_equal(fitted(mySBM), predict(mySBM))

  ## correctness
  expect_lt(rmse(sort(mySBM$connectParam$mean), sort(means)), 1e-1)
  expect_lt(1 - aricode::ARI(mySBM$memberships[[1]], mySampler$memberships[[1]]), 1e-1)
  expect_lt(1 - aricode::ARI(mySBM$memberships[[2]], mySampler$memberships[[2]]), 1e-1)

  ## prediction wrt BM
  for (Q in mySBM$storedModels$indexModel) {
    pred_bm  <- BM_out$prediction(Q = Q)
    mySBM$setModel(Q-1)
    pred_sbm <- predict(mySBM)
    expect_lt( rmse(pred_bm, pred_sbm), 1e-12)
  }

})

test_that("active bindings are working in the class", {
  A <- matrix(rbinom(200,1,.2),20,10)

  myBipartite <- BipartiteSBM_fit$new(incidenceMatrix = A,model = "bernoulli",dimLabels = c("Actor","Stuff"))

  tau1 <- matrix(runif(20*2),20,2)
  tau1 <- tau1 / rowSums(tau1)
  tau2 <- matrix(runif(10*3),10,3)
  tau2 <- tau2 / rowSums(tau2)
  myBipartite$probMemberships <- list(tau1,tau2)
  myBipartite$blockProp <- list(colMeans(tau1),colMeans(tau2))
  myBipartite$connectParam <- list(mean = matrix(runif(3*2),3,2))

  expect_equal(unname(myBipartite$nbNodes),c(20,10))

  expect_equal(myBipartite$memberships[[1]], 1+(tau1[,1]<.5)*1)
  expect_equal(dim(myBipartite$connectParam$mean),c(3,2))
  expect_equal(length(myBipartite$blockProp),2)

})
