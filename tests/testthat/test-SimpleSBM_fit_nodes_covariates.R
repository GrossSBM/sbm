set.seed(1234)

rmse <- function(theta, theta_star) {
    sqrt(sum((theta - theta_star)^2) / sum(theta_star^2))
}

## Common parameters
nbNodes <- 40
nbBlocks <- 3
blockProp <- c(1 / 3, 1 / 3, 1 / 3)

# Nodes covariates parameters
nbNodesCovar <- 2
nodesCovar <- matrix(rnorm(nbNodes * nbNodesCovar), nbNodes, nbNodesCovar)
nodesCovarParam <- matrix(rnorm(nbBlocks * nbNodesCovar), nbBlocks, nbNodesCovar)

# Edge covariates parameters
covarParam <- c(-2, 2)
covar1 <- matrix(rnorm(nbNodes**2), nbNodes, nbNodes)
covar2 <- matrix(rnorm(nbNodes**2), nbNodes, nbNodes)
covar1 <- covar1 + t(covar1)
covar2 <- covar2 + t(covar2)
covarList <- list(covar1 = covar1, covar2 = covar2)

test_that("SimpleSBM_fit 'Bernoulli' model, undirected, with nodes covariates", {
    ## SIMPLE UNDIRECTED BERNOULLI SBM WITH NODES COVARIATES
    means <- diag(.4, 3) + 0.05
    connectParam <- list(mean = means)

    ## Basic construction with nodes covariates
    mySampler <- SimpleSBM$new("bernoulli", nbNodes, FALSE, blockProp, connectParam,
        nodesCovar = nodesCovar
    )
    mySampler$rMemberships(store = TRUE)
    mySampler$rEdges(store = TRUE)

    ## Construction----------------------------------------------------------------
    mySBM <- SimpleSBM_fit$new(mySampler$networkData, "bernoulli", FALSE,
        nodesCovar = nodesCovar
    )

    ## Checking class
    expect_true(inherits(mySBM, "SBM"))
    expect_true(inherits(mySBM, "SimpleSBM"))
    expect_true(inherits(mySBM, "SimpleSBM_fit"))

    ## Checking field access prior to estimation
    ## parameters
    expect_equal(mySBM$modelName, "bernoulli")
    expect_equal(unname(mySBM$nbNodes), nbNodes)
    expect_equal(mySBM$dimLabels, c(node = "nodeName"))

    ## nodes covariates
    expect_equal(mySBM$nbNodesCovariates, c(node = nbNodesCovar))
    expect_equal(dim(mySBM$nodesCovariates[["node"]]), c(nbNodes, nbNodesCovar))
    expect_equal(mySBM$nodesCovariates, list(node = nodesCovar))

    ## Estimation-----------------------------------------------------------------
    BM_out <- mySBM$optimize(estimOptions = list(verbosity = 0, fast = TRUE))
    mySBM$setModel(nbBlocks)

    ## Expectation
    expect_equal(dim(mySBM$expectation), c(nbNodes, nbNodes))
    expect_true(all(mySBM$expectation >= 0, na.rm = TRUE))
    expect_true(all(mySBM$expectation <= 1, na.rm = TRUE))

    ## blocks
    expect_equal(mySBM$nbBlocks, nbBlocks)
    expect_equal(dim(mySBM$probMemberships), c(nbNodes, nbBlocks))
    expect_equal(sort(unique(mySBM$memberships)), 1:nbBlocks)

    ## nodes covariate parameters should be initialized
    expect_true(!is.null(mySBM$nodesCovarParam))
    expect_true(length(mySBM$nodesCovarParam) > 0)
})

test_that("SimpleSBM_fit 'Bernoulli' model, directed, with nodes covariates", {
    ## SIMPLE DIRECTED BERNOULLI SBM WITH NODES COVARIATES
    means <- matrix(c(
        0.45, 0.10, 0.05,
        0.20, 0.40, 0.10,
        0.15, 0.05, 0.35
    ), 3, 3, byrow = TRUE)
    connectParam <- list(mean = means)

    ## Basic construction
    mySampler <- SimpleSBM$new("bernoulli", nbNodes, TRUE, blockProp, connectParam,
        nodesCovar = nodesCovar
    )
    mySampler$rMemberships(store = TRUE)
    mySampler$rEdges(store = TRUE)

    ## Construction
    mySBM <- SimpleSBM_fit$new(mySampler$networkData, "bernoulli", TRUE,
        nodesCovar = nodesCovar
    )

    ## Checking class
    expect_true(inherits(mySBM, "SBM"))
    expect_true(inherits(mySBM, "SimpleSBM_fit"))

    ## Checking field access prior to estimation
    expect_equal(mySBM$directed, TRUE)
    expect_equal(mySBM$nbNodesCovariates, c(node = nbNodesCovar))
    expect_equal(mySBM$nbDyads, nbNodes * (nbNodes - 1))

    ## Estimation
    BM_out <- mySBM$optimize(estimOptions = list(verbosity = 0, fast = TRUE))
    mySBM$setModel(nbBlocks)

    ## Expectation
    expect_equal(dim(mySBM$expectation), c(nbNodes, nbNodes))
    expect_true(all(mySBM$expectation >= 0, na.rm = TRUE))
    expect_true(all(mySBM$expectation <= 1, na.rm = TRUE))
})

test_that("SimpleSBM_fit 'Poisson' model, undirected, with nodes covariates", {
    ## SIMPLE UNDIRECTED POISSON SBM WITH NODES COVARIATES
    means <- matrix(c(
        12, 4, 3,
        4, 11, 2,
        3, 2, 10
    ), 3, 3, byrow = TRUE)
    connectParam <- list(mean = means)

    ## Basic construction
    mySampler <- SimpleSBM$new("poisson", nbNodes, FALSE, blockProp, connectParam,
        nodesCovar = nodesCovar
    )
    mySampler$rMemberships(store = TRUE)
    mySampler$rEdges(store = TRUE)

    ## Construction
    mySBM <- SimpleSBM_fit$new(mySampler$networkData, "poisson", FALSE,
        nodesCovar = nodesCovar
    )

    ## Checking class
    expect_true(inherits(mySBM, "SimpleSBM_fit"))

    ## Checking nodes covariates
    expect_equal(mySBM$nbNodesCovariates, c(node = nbNodesCovar))
    expect_equal(nrow(mySBM$nodesCovariates[["node"]]), nbNodes)
    expect_equal(ncol(mySBM$nodesCovariates[["node"]]), nbNodesCovar)

    ## Estimation
    BM_out <- mySBM$optimize(estimOptions = list(verbosity = 0, fast = TRUE))
    mySBM$setModel(nbBlocks)

    ## Expectation
    expect_equal(dim(mySBM$expectation), c(nbNodes, nbNodes))
    expect_true(all(mySBM$expectation >= 0, na.rm = TRUE))
})

test_that("SimpleSBM_fit without nodes covariates returns zero dimensions", {
    ## SIMPLE UNDIRECTED BERNOULLI SBM WITHOUT NODES COVARIATES
    means <- diag(.4, 3) + 0.05
    connectParam <- list(mean = means)

    ## Basic construction without nodes covariates
    mySampler <- SimpleSBM$new("bernoulli", nbNodes, FALSE, blockProp, connectParam)
    mySampler$rMemberships(store = TRUE)
    mySampler$rEdges(store = TRUE)

    ## Construction without nodes covariates
    mySBM <- SimpleSBM_fit$new(mySampler$networkData, "bernoulli", FALSE)

    ## Check that nbNodesCovariates returns named vector of zeros
    expect_equal(mySBM$nbNodesCovariates, c(node = 0))

    ## nodesCovariates should be empty
    expect_equal(length(mySBM$nodesCovariates), 0)
})
