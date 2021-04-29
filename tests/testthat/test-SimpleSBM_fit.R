set.seed(1234)

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
  mySampler <- SimpleSBM$new('bernoulli', nbNodes, FALSE, blockProp, connectParam)
  mySampler$rMemberships(store = TRUE)
  mySampler$rEdges(store = TRUE)

  ## Construction----------------------------------------------------------------
  mySBM <- SimpleSBM_fit$new(mySampler$networkData, 'bernoulli', FALSE)
  expect_error(SimpleSBM_fit$new(mySampler$networkData, 'bernouilli', FALSE))
  expect_error(SimpleSBM_fit$new(mySampler$networkData[1:20, 1:30], 'bernouilli', FALSE))
  expect_error(SimpleSBM_fit$new(mySampler$networkData, 'bernoulli', TRUE))

  ## Checking class
  expect_true(inherits(mySBM, "SBM"))
  expect_true(inherits(mySBM, "SimpleSBM"))
  expect_true(inherits(mySBM, "SimpleSBM_fit"))

  ## Checking field access and format prior to estimation
  ## parameters
  expect_equal(mySBM$modelName, 'bernoulli')
  expect_equal(unname(mySBM$nbNodes), nbNodes)
  expect_equal(mySBM$nbDyads, nbNodes*(nbNodes - 1)/2)
  expect_true(all(is.na(diag(mySBM$networkData))))
  expect_true(isSymmetric(mySBM$networkData))
  expect_true(!mySBM$directed)
  expect_true(is.matrix(mySBM$connectParam$mean))

  ## covariates
  expect_null(mySBM$covarExpect)
  expect_equal(mySBM$nbCovariates, 0)
  expect_equal(mySBM$covarList, list())
  expect_equal(mySBM$covarParam, numeric(0))
  expect_equal(mySBM$covarEffect, numeric(0))

  ## S3 methods
  expect_equal(coef(mySBM, 'connectivity'), mySBM$connectParam)
  expect_equal(coef(mySBM, 'block')       , mySBM$blockProp)
  expect_equal(coef(mySBM, 'covariates')  , mySBM$covarParam)

  ## Estimation-----------------------------------------------------------------
  BM_out <- mySBM$optimize(estimOptions=list(verbosity = 0))
  mySBM$setModel(3)

  ## Field set after optimization
  expect_equal(mySBM$nbConnectParam, nbBlocks * (nbBlocks + 1)/2 )
  expect_equal(mySBM$penalty, (nbBlocks * (nbBlocks + 1))/2 * log(nbNodes *(nbNodes-1)/2) +  (nbBlocks - 1) * log(nbNodes))
  expect_equal(mySBM$entropy, -sum(mySBM$probMemberships * log(mySBM$probMemberships)))

  ## Expectation
  expect_equal(dim(mySBM$expectation), c(nbNodes, nbNodes))
  expect_true(all(mySBM$expectation >= 0, na.rm = TRUE))
  expect_true(all(mySBM$expectation <= 1, na.rm = TRUE))
  expect_null(mySBM$connectParam$var)

  ## S3 methods
  expect_equal(coef(mySBM, 'connectivity'), mySBM$connectParam)
  expect_equal(coef(mySBM, 'block')       , mySBM$blockProp)
  expect_equal(coef(mySBM, 'covariates')  , mySBM$covarParam)
  expect_equal(mySBM$predict(), predict(mySBM))
  expect_equal(fitted(mySBM), predict(mySBM))

  ## blocks
  expect_equal(mySBM$nbBlocks, nbBlocks)
  expect_equal(dim(mySBM$probMemberships), c(nbNodes, nbBlocks))
  expect_equal(sort(unique(mySBM$memberships)), 1:nbBlocks)
  expect_equal(length(mySBM$memberships), nbNodes)

  ## correctness
  expect_lt(rmse(mySBM$connectParam$mean, means), 0.25)
  expect_lt(1 - aricode::ARI(mySBM$memberships, mySampler$memberships), 0.25)

  ## prediction wrt BM
  for (Q in mySBM$storedModels$indexModel) {
    pred_bm  <- BM_out$prediction(Q = Q)
    mySBM$setModel(Q)
    pred_sbm <- predict(mySBM)
    expect_lt( rmse(pred_bm, pred_sbm), 1e-12)
  }

})

test_that("SimpleSBM_fit 'Bernoulli' model, directed, no covariate", {

  ## SIMPLE UNDIRECTED BERNOULLI SBM
  means <- matrix(c(1:9)/10, 3,  3)
  connectParam <- list(mean = means)

  ## Basic construction - check for wrong specifications
  mySampler <- SimpleSBM$new('bernoulli', nbNodes, TRUE, blockProp, connectParam)
  mySampler$rMemberships(store = TRUE)
  mySampler$rEdges(store = TRUE)

  ## Construction----------------------------------------------------------------
  mySBM <- SimpleSBM_fit$new(mySampler$networkData, 'bernoulli', TRUE)
  expect_error(SimpleSBM_fit$new(mySampler$networkData, 'bernouilli', TRUE))
  expect_error(SimpleSBM_fit$new(mySampler$networkData, 'bernoulli', FALSE))
  expect_error(SimpleSBM_fit$new(mySampler$networkData[1:20, 1:30], 'bernouilli', FALSE))

  ## Checking class
  expect_true(inherits(mySBM, "SBM"))
  expect_true(inherits(mySBM, "SimpleSBM"))
  expect_true(inherits(mySBM, "SimpleSBM_fit"))

  ## Checking field access and format prior to estimation
  ## parameters
  expect_equal(mySBM$modelName, 'bernoulli')
  expect_equal(unname(mySBM$nbNodes), nbNodes)
  expect_equal(mySBM$nbDyads, nbNodes*(nbNodes - 1))
  expect_true(all(is.na(diag(mySBM$networkData))))
  expect_true(!isSymmetric(mySBM$networkData))
  expect_true(mySBM$directed)
  expect_true(is.matrix(mySBM$connectParam$mean))

  ## covariates
  expect_null(mySBM$covarExpect)
  expect_equal(mySBM$nbCovariates, 0)
  expect_equal(mySBM$covarList, list())
  expect_equal(mySBM$covarParam, numeric(0))
  expect_equal(mySBM$covarEffect, numeric(0))

  ## S3 methods
  expect_equal(coef(mySBM, 'connectivity'), mySBM$connectParam)
  expect_equal(coef(mySBM, 'block')       , mySBM$blockProp)
  expect_equal(coef(mySBM, 'covariates')  , mySBM$covarParam)

  ## Estimation-----------------------------------------------------------------
  BM_out <- mySBM$optimize(estimOptions=list(verbosity = 0))
  mySBM$setModel(3)

  expect_equal(mySBM$nbConnectParam, nbBlocks * nbBlocks)
  expect_equal(mySBM$penalty, nbBlocks * nbBlocks * log(nbNodes * (nbNodes -1)) +  (nbBlocks - 1) * log(nbNodes))
  expect_equal(mySBM$entropy, -sum(mySBM$probMemberships * log(mySBM$probMemberships)))

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
  expect_equal(fitted(mySBM), predict(mySBM))

  ## correctness
  expect_lt(rmse(sort(mySBM$connectParam$mean), means), 0.2)
  expect_lt(1 - aricode::ARI(mySBM$memberships, mySampler$memberships), 0.2)

  ## prediction wrt BM
  for (Q in mySBM$storedModels$indexModel) {
    pred_bm  <- BM_out$prediction(Q = Q)
    mySBM$setModel(Q)
    pred_sbm <- predict(mySBM)
    expect_lt( rmse(pred_bm, pred_sbm), 1e-12)
  }

})

test_that("SimpleSBM_fit 'Poisson' model, undirected, no covariate", {

  ## SIMPLE UNDIRECTED POISSON SBM
  means <- diag(15., 3) + 5
  connectParam <- list(mean = means)

  ## Basic construction - check for wrong specifications
  mySampler <- SimpleSBM$new('poisson', nbNodes, FALSE, blockProp, connectParam)
  mySampler$rMemberships(store = TRUE)
  mySampler$rEdges(store = TRUE)

  ## Construction----------------------------------------------------------------
  mySBM <- SimpleSBM_fit$new(mySampler$networkData, 'poisson', FALSE)
  expect_error(SimpleSBM_fit$new(SamplerBernoulli$networkData, 'poison', FALSE))
  expect_error(SimpleSBM_fit$new(SamplerBernoulli$networkData[1:20, 1:30], 'poisson', FALSE))
  expect_error(SimpleSBM_fit$new(SamplerBernoulli$networkData, 'poisson', TRUE))

  ## Checking class
  expect_true(inherits(mySBM, "SBM"))
  expect_true(inherits(mySBM, "SimpleSBM"))
  expect_true(inherits(mySBM, "SimpleSBM_fit"))

  ## Checking field access and format prior to estimation
  ## parameters
  expect_equal(mySBM$modelName, 'poisson')
  expect_equal(unname(mySBM$nbNodes), nbNodes)
  expect_equal(mySBM$nbDyads, nbNodes*(nbNodes - 1)/2)
  expect_true(all(is.na(diag(mySBM$networkData))))
  expect_true(isSymmetric(mySBM$networkData))
  expect_true(!mySBM$directed)
  expect_true(is.matrix(mySBM$connectParam$mean))

  ## covariates
  expect_null(mySBM$covarExpect)
  expect_equal(mySBM$nbCovariates, 0)
  expect_equal(mySBM$covarList, list())
  expect_equal(mySBM$covarParam, numeric(0))
  expect_equal(mySBM$covarEffect, numeric(0))

  ## S3 methods
  expect_equal(coef(mySBM, 'connectivity'), mySBM$connectParam)
  expect_equal(coef(mySBM, 'block')       , mySBM$blockProp)
  expect_equal(coef(mySBM, 'covariates')  , mySBM$covarParam)

  ## Estimation-----------------------------------------------------------------
  BM_out <- mySBM$optimize(estimOptions=list(verbosity = 0))
  mySBM$setModel(3)

  ## Expectation
  expect_equal(dim(mySBM$expectation), c(nbNodes, nbNodes))
  expect_true(all(mySampler$expectation >= 0, na.rm = TRUE))
  expect_null(mySBM$connectParam$var)

  ## S3 methods
  expect_equal(coef(mySBM, 'connectivity'), mySBM$connectParam)
  expect_equal(coef(mySBM, 'block')       , mySBM$blockProp)
  expect_equal(coef(mySBM, 'covariates')  , mySBM$covarParam)
  expect_equal(mySBM$predict(), predict(mySBM))
  expect_equal(fitted(mySBM), predict(mySBM))

  ## blocks
  expect_equal(mySBM$nbBlocks, nbBlocks)
  expect_equal(dim(mySBM$probMemberships), c(nbNodes, nbBlocks))
  expect_equal(sort(unique(mySBM$memberships)), 1:nbBlocks)
  expect_equal(length(mySBM$memberships), nbNodes)

  ## correctness
  expect_lt(rmse(mySBM$connectParam$mean, means), 1e-1)
  expect_lt(1 - aricode::ARI(mySBM$memberships, mySampler$memberships), 1e-1)

  ## prediction wrt BM
  for (Q in mySBM$storedModels$indexModel) {
    pred_bm  <- BM_out$prediction(Q = Q)
    mySBM$setModel(Q)
    pred_sbm <- predict(mySBM)
    expect_lt( rmse(pred_bm, pred_sbm), 1e-12)
  }

})

test_that("SimpleSBM_fit 'Poisson' model, directed, no covariate", {

  ## SIMPLE UNDIRECTED BERNOULLI SBM
  means <- matrix(1:9, 3,  3)
  connectParam <- list(mean = means)

  ## Basic construction - check for wrong specifications
  mySampler <- SimpleSBM$new('poisson', nbNodes, TRUE, blockProp, connectParam)
  mySampler$rMemberships(store = TRUE)
  mySampler$rEdges(store = TRUE)

  ## Construction----------------------------------------------------------------
  mySBM <- SimpleSBM_fit$new(mySampler$networkData, 'poisson', TRUE)
  expect_error(SimpleSBM_fit$new(SamplerBernoulli$networkData, 'poison', TRUE))
  expect_error(SimpleSBM_fit$new(SamplerBernoulli$networkData, 'poisson', FALSE))

  ## Checking class
  expect_true(inherits(mySBM, "SBM"))
  expect_true(inherits(mySBM, "SimpleSBM"))
  expect_true(inherits(mySBM, "SimpleSBM_fit"))

  ## Checking field access and format prior to estimation
  ## parameters
  expect_equal(mySBM$modelName, 'poisson')
  expect_equal(unname(mySBM$nbNodes), nbNodes)
  expect_equal(mySBM$nbDyads, nbNodes*(nbNodes - 1))
  expect_true(all(is.na(diag(mySBM$networkData))))
  expect_true(!isSymmetric(mySBM$networkData))
  expect_true(mySBM$directed)
  expect_true(is.matrix(mySBM$connectParam$mean))

  ## covariates
  expect_null(mySBM$covarExpect)
  expect_equal(mySBM$nbCovariates, 0)
  expect_equal(mySBM$covarList, list())
  expect_equal(mySBM$covarParam, numeric(0))
  expect_equal(mySBM$covarEffect, numeric(0))

  ## S3 methods
  expect_equal(coef(mySBM, 'connectivity'), mySBM$connectParam)
  expect_equal(coef(mySBM, 'block')       , mySBM$blockProp)
  expect_equal(coef(mySBM, 'covariates')  , mySBM$covarParam)

  ## Estimation-----------------------------------------------------------------
  BM_out <- mySBM$optimize(estimOptions=list(verbosity = 0))
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

  ## S3 methods
  expect_equal(coef(mySBM, 'connectivity'), mySBM$connectParam)
  expect_equal(coef(mySBM, 'block')       , mySBM$blockProp)
  expect_equal(coef(mySBM, 'covariates')  , mySBM$covarParam)
  expect_equal(mySBM$predict(), predict(mySBM))
  expect_equal(fitted(mySBM), predict(mySBM))

  ## correctness
  expect_lt(rmse(sort(mySBM$connectParam$mean), means), 1e-1)
  expect_lt(1 - aricode::ARI(mySBM$memberships, mySampler$memberships), 1e-1)

  ## prediction wrt BM
  for (Q in mySBM$storedModels$indexModel) {
    pred_bm  <- BM_out$prediction(Q = Q)
    mySBM$setModel(Q)
    pred_sbm <- predict(mySBM)
    expect_lt( rmse(pred_bm, pred_sbm), 1e-12)
  }

})


test_that("SimpleSBM_fit 'Gaussian' model, undirected, no covariate", {

  ## SIMPLE UNDIRECTED GAUSSIAN SBM
  means <- diag(15., 3) + 5 # connectivity matrix: affiliation network
  connectParam <- list(mean = means, var = 2)

  ## Basic construction - check for wrong specifications
  mySampler <- SimpleSBM$new('gaussian', nbNodes, FALSE, blockProp, connectParam)
  mySampler$rMemberships(store = TRUE)
  mySampler$rEdges(store = TRUE)

  ## Construction----------------------------------------------------------------
  mySBM <- SimpleSBM_fit$new(mySampler$networkData, 'gaussian', FALSE)
  expect_error(SimpleSBM_fit$new(SamplerBernoulli$networkData, 'normal', FALSE))
  expect_error(SimpleSBM_fit$new(SamplerBernoulli$networkData[1:20, 1:30], 'gaussian', FALSE))
  expect_error(SimpleSBM_fit$new(SamplerBernoulli$networkData, 'gaussian', TRUE))

  ## Checking class
  expect_true(inherits(mySBM, "SBM"))
  expect_true(inherits(mySBM, "SimpleSBM"))
  expect_true(inherits(mySBM, "SimpleSBM_fit"))

  ## Checking field access and format prior to estimation
  ## parameters
  expect_equal(mySBM$modelName, 'gaussian')
  expect_equal(unname(mySBM$nbNodes), nbNodes)
  expect_equal(mySBM$nbDyads, nbNodes*(nbNodes - 1)/2)
  expect_true(all(is.na(diag(mySBM$networkData))))
  expect_true(isSymmetric(mySBM$networkData))
  expect_true(!mySBM$directed)
  expect_true(is.matrix(mySBM$connectParam$mean))

  ## covariates
  expect_null(mySBM$covarExpect)
  expect_equal(mySBM$nbCovariates, 0)
  expect_equal(mySBM$covarList, list())
  expect_equal(mySBM$covarParam, numeric(0))

  ## S3 methods
  expect_equal(coef(mySBM, 'connectivity'), mySBM$connectParam)
  expect_equal(coef(mySBM, 'block')       , mySBM$blockProp)
  expect_equal(coef(mySBM, 'covariates')  , mySBM$covarParam)

  ## Estimation-----------------------------------------------------------------
  BM_out <- mySBM$optimize(estimOptions=list(verbosity = 0))
  mySBM$setModel(3)

  ## Expectation
  expect_equal(dim(mySBM$expectation), c(nbNodes, nbNodes))
  expect_gt(mySBM$connectParam$var, 0)

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
  expect_equal(fitted(mySBM), predict(mySBM))

  ## correctness
  expect_lt(rmse(mySBM$connectParam$mean, means), 1e-1)
  expect_lt(1 - aricode::ARI(mySBM$memberships, mySampler$memberships), 1e-1)

  ## prediction wrt BM
  for (Q in mySBM$storedModels$indexModel) {
    pred_bm  <- BM_out$prediction(Q = Q)
    mySBM$setModel(Q)
    pred_sbm <- predict(mySBM)
    expect_lt( rmse(pred_bm, pred_sbm), 1e-12)
  }

})

test_that("SimpleSBM_fit 'Gaussian' model, undirected, no covariate", {

  ## SIMPLE UNDIRECTED GAUSSIAN SBM
  means <- matrix(1:9,3,3)
  connectParam <- list(mean = means, var = 2)

  ## Basic construction - check for wrong specifications
  mySampler <- SimpleSBM$new('gaussian', nbNodes, TRUE, blockProp, connectParam)
  mySampler$rMemberships(store = TRUE)
  mySampler$rEdges(store = TRUE)

  ## Construction----------------------------------------------------------------
  mySBM <- SimpleSBM_fit$new(mySampler$networkData, 'gaussian', TRUE)
  expect_error(SimpleSBM_fit$new(SamplerBernoulli$networkData, 'normal', TRUE))
  expect_error(SimpleSBM_fit$new(SamplerBernoulli$networkData[1:20, 1:30], 'gaussian', TRUE))
  expect_error(SimpleSBM_fit$new(SamplerBernoulli$networkData, 'gaussian', FALSE))

  ## Checking class
  expect_true(inherits(mySBM, "SBM"))
  expect_true(inherits(mySBM, "SimpleSBM"))
  expect_true(inherits(mySBM, "SimpleSBM_fit"))

  ## Checking field access and format prior to estimation
  ## parameters
  expect_equal(mySBM$modelName, 'gaussian')
  expect_equal(unname(mySBM$nbNodes), nbNodes)
  expect_equal(mySBM$nbDyads, nbNodes*(nbNodes - 1))
  expect_true(all(is.na(diag(mySBM$networkData))))
  expect_true(!isSymmetric(mySBM$networkData))
  expect_true(mySBM$directed)
  expect_true(is.matrix(mySBM$connectParam$mean))

  ## covariates
  expect_null(mySBM$covarExpect)
  expect_equal(mySBM$nbCovariates, 0)
  expect_equal(mySBM$covarList, list())
  expect_equal(mySBM$covarParam, numeric(0))
  expect_equal(mySBM$covarEffect, numeric(0))

  ## S3 methods
  expect_equal(coef(mySBM, 'connectivity'), mySBM$connectParam)
  expect_equal(coef(mySBM, 'block')       , mySBM$blockProp)
  expect_equal(coef(mySBM, 'covariates')  , mySBM$covarParam)

  ## Estimation-----------------------------------------------------------------
  BM_out <- mySBM$optimize(estimOptions=list(verbosity = 0))
  mySBM$setModel(3)

  ## Expectation
  expect_equal(dim(mySBM$expectation), c(nbNodes, nbNodes))
  expect_gt(mySBM$connectParam$var, 0)

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
  expect_equal(fitted(mySBM), predict(mySBM))

  ## correctness
  expect_lt(rmse(sort(mySBM$connectParam$mean), means), 1e-1)
  expect_lt(1 - aricode::ARI(mySBM$memberships, mySampler$memberships), 1e-1)

  ## prediction wrt BM
  for (Q in mySBM$storedModels$indexModel) {
    pred_bm  <- BM_out$prediction(Q = Q)
    mySBM$setModel(Q)
    pred_sbm <- predict(mySBM)
    expect_lt( rmse(pred_bm, pred_sbm), 1e-12)
  }

})

test_that("active binding are working in the class", {

  A <- matrix(rbinom(100,1,.2),10,10)

  mySimple <- SimpleSBM_fit$new(adjacencyMatrix = A,model = "bernoulli",directed = TRUE, dimLabels = "Actor")
  p <- runif(10)
  mySimple$probMemberships <- matrix(c(p,1-p),10,2,byrow=FALSE)
  mySimple$blockProp <- c(.3,.7)
  mySimple$connectParam <- list(mean = matrix(runif(4),2,2))
  expect_equal(mySimple$memberships, 1+(p<.5)*1)
  expect_equal(dim(mySimple$connectParam$mean),c(2,2))
  expect_equal(length(mySimple$blockProp),2)

})

