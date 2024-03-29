if (Sys.info()['sysname'] != "Windows") {
  set.seed(123)

  rmse <- function(theta, theta_star) { sqrt(sum((theta - theta_star)^2)/sum(theta_star^2)) }

  ## Common parameters
  nbNodes  <- c(30, 60)
  blockProp <- list(row = c(.5, .5), col = c(1/3, 1/3, 1/3)) # group proportions
  nbBlocks <- sapply(blockProp, length)
  covarParam <- c(-2,2)
  covar1 <- matrix(rnorm(prod(nbNodes)), nbNodes[1], nbNodes[2])
  covar2 <- matrix(rnorm(prod(nbNodes)), nbNodes[1], nbNodes[2])
  covarList <- list(covar1 = covar1, covar2 = covar2)

  test_that("BipartiteSBM_fit 'Bernoulli' model, undirected, no covariate", {

    ## BIPARTITE UNDIRECTED BERNOULLI SBM
    means <- matrix(c(0.05, 0.95, 0.4, 0.75, 0.15, 0.6), 2, 3)  # connectivity matrix
    connectParam <- list(mean = means)

    ## Basic construction - check for wrong specifications
    mySampler <- BipartiteSBM$new('bernoulli', nbNodes, blockProp, connectParam,covarParam = covarParam[1], covarList = covarList[1])
    mySampler$rMemberships(store = TRUE)
    mySampler$rEdges(store = TRUE)

    ## Construction----------------------------------------------------------------
    mySBM <- BipartiteSBM_fit$new(mySampler$networkData, 'bernoulli', covarList = covarList[1])
    expect_error(BipartiteSBM_fit$new(SamplerBernoulli$networkData, 'bernouilli', covarList = covarList[1]))

    ## Checking class
    expect_true(inherits(mySBM, "SBM"))
    expect_true(inherits(mySBM, "BipartiteSBM"))
    expect_true(inherits(mySBM, "BipartiteSBM_fit"))

    ## Checking field access and format prior to estimation
    ## parameters
    expect_equal(mySBM$modelName, 'bernoulli')
    expect_equal(unname(mySBM$nbNodes), nbNodes)
    expect_equal(mySBM$dimLabels, c(row="row", col="col"))
    expect_equal(mySBM$nbDyads, nbNodes[1]*nbNodes[2])
    expect_true(is.matrix(mySBM$connectParam$mean))

    ## covariates
    expect_equal(dim(mySBM$covarEffect), nbNodes)
    expect_equal(mySBM$nbCovariates, 1)
    expect_equal(mySBM$covarList, covarList[1])
    expect_equal(mySBM$covarParam, 0)

    ## S3 methods
    expect_equal(coef(mySBM, 'connectivity'), mySBM$connectParam)
    expect_equal(coef(mySBM, 'block')       , mySBM$blockProp)
    expect_equal(coef(mySBM, 'covariates')  , mySBM$covarParam)

    ## Estimation-----------------------------------------------------------------
    BM_out <- mySBM$optimize(estimOptions  = list(verbosity = 0, fast = TRUE))
    mySBM$setModel(4)

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
    expect_equal(predict(mySBM, covarList[1]), fitted(mySBM))
    expect_error(predict(mySBM, covarList))

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
    mySampler <- BipartiteSBM$new('poisson', nbNodes, blockProp, connectParam, covarParam = covarParam, covarList = covarList)
    mySampler$rMemberships(store = TRUE)
    mySampler$rEdges(store = TRUE)

    ## Construction----------------------------------------------------------------
    mySBM <- BipartiteSBM_fit$new(mySampler$networkData, 'poisson', covarList = covarList)
    expect_error(BipartiteSBM_fit$new(SamplerBernoulli$networkData, 'poison', covarList = covarList))

    ## Checking class
    expect_true(inherits(mySBM, "SBM"))
    expect_true(inherits(mySBM, "BipartiteSBM"))
    expect_true(inherits(mySBM, "BipartiteSBM_fit"))

    ## Checking field access and format prior to estimation
    ## parameters
    expect_equal(mySBM$modelName, 'poisson')
    expect_equal(unname(mySBM$nbNodes), nbNodes)
    expect_equal(mySBM$dimLabels, c(row="row", col="col"))
    expect_equal(mySBM$nbDyads, nbNodes[1]*nbNodes[2])
    expect_true(is.matrix(mySBM$connectParam$mean))

    ## covariates
    expect_equal(dim(mySBM$covarEffect), nbNodes)
    expect_equal(mySBM$nbCovariates, 2)
    expect_equal(mySBM$covarList, covarList)
    expect_equal(mySBM$covarParam, c(0,0))

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

    ## correctness
    expect_lt(rmse(sort(mySBM$connectParam$mean), sort(means)), 1e-1)
    expect_lt(1 - aricode::ARI(mySBM$memberships[[1]], mySampler$memberships[[1]]), 1e-1)
    expect_lt(1 - aricode::ARI(mySBM$memberships[[2]], mySampler$memberships[[2]]), 1e-1)

    ## S3 methods
    expect_equal(coef(mySBM, 'connectivity'), mySBM$connectParam)
    expect_equal(coef(mySBM, 'block')       , mySBM$blockProp)
    expect_equal(coef(mySBM, 'covariates')  , mySBM$covarParam)
    expect_equal(mySBM$predict(), predict(mySBM))
    expect_equal(predict(mySBM, covarList), fitted(mySBM))
    expect_error(predict(mySBM, covarList[1]))

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
    mySampler <- BipartiteSBM$new('gaussian', nbNodes, blockProp, connectParam, covarParam = covarParam, covarList = covarList)
    mySampler$rMemberships(store = TRUE)
    mySampler$rEdges(store = TRUE)

    ## Construction----------------------------------------------------------------
    mySBM <- BipartiteSBM_fit$new(mySampler$networkData, 'gaussian', covarList = covarList)
    expect_error(BipartiteSBM_fit$new(SamplerBernoulli$networkData, 'groß', covarList = covarList))

    ## Checking class
    expect_true(inherits(mySBM, "SBM"))
    expect_true(inherits(mySBM, "BipartiteSBM"))
    expect_true(inherits(mySBM, "BipartiteSBM_fit"))

    ## Checking field access and format prior to estimation
    ## parameters
    expect_equal(mySBM$modelName, 'gaussian')
    expect_equal(unname(mySBM$nbNodes), nbNodes)
    expect_equal(mySBM$dimLabels, c(row="row", col="col"))
    expect_equal(mySBM$nbDyads, nbNodes[1]*nbNodes[2])
    expect_true(is.matrix(mySBM$connectParam$mean))

    ## covariates
    expect_equal(dim(mySBM$covarEffect), nbNodes)
    expect_equal(mySBM$nbCovariates, 2)
    expect_equal(mySBM$covarList, covarList)
    expect_equal(mySBM$covarParam, c(0,0))

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
    expect_equal(predict(mySBM, covarList), fitted(mySBM))

    ## correctness
    expect_lt(rmse(sort(mySBM$connectParam$mean), sort(means)), 1e-1)
    expect_lt(1 - aricode::ARI(mySBM$memberships[[1]], mySampler$memberships[[1]]), 2e-1)
    expect_lt(1 - aricode::ARI(mySBM$memberships[[2]], mySampler$memberships[[2]]), 2e-1)

    ## prediction wrt BM
    for (Q in mySBM$storedModels$indexModel) {
      pred_bm  <- BM_out$prediction(Q = Q)
      mySBM$setModel(Q-1)
      pred_sbm <- predict(mySBM)
      expect_lt( rmse(pred_bm, pred_sbm), 1e-12)
    }

  })

}
