rm(list=ls())
library(sbm)
set.seed(123)

rmse <- function(theta, theta_star) {
    sqrt(sum((theta - theta_star)^2) / sum(theta_star^2))
}

## Common parameters
nbNodes <- c(30, 60)
blockProp <- list(row = c(1 / 3, 1 / 3, 1 / 3), col = c(1 / 3, 1 / 3, 1 / 3))
nbBlocks <- sapply(blockProp, length)

# Nodes covariates parameters
nbNodesCovarRow <- 2+1
nbNodesCovarCol <- 3+1
nodesCovarRow <- cbind(1,matrix(rnorm(nbNodes[1] * (nbNodesCovarRow-1)), nbNodes[1], nbNodesCovarRow-1))
nodesCovarCol <- cbind(1,matrix(rnorm(nbNodes[2] * (nbNodesCovarCol-1)), nbNodes[2], nbNodesCovarCol-1))
nodesCovarList <- list(row = nodesCovarRow, col = nodesCovarCol)


B <- matrix(0,nbNodesCovarRow,nbBlocks[1])
B[1,] <- c(-0.01, 3, 4)
B[2,] <- c( 2,    0.4,  -2 )
B[3,] <- c(1, 1, 3)
B_contrQ <- B - B[,nbBlocks[1]]

G <- matrix(0,nbNodesCovarCol,nbBlocks[2])
G[1,] <- c(-0, 1, 0)
G[2,] <- c( 2, 1, 0.4)
G[3,] <- c( -2, 1, 2)
G[4,] <- c( 0, 0, 2)
G_contrQ <- G - G[,nbBlocks[2]]

nodesCovariatesParam = vector('list',2)
nodesCovariatesParam[[1]] <- B_contrQ
nodesCovariatesParam[[2]] <- G_contrQ


# # Edge covariates parameters
# covarParam <- c(-2, 2)
# covar1 <- matrix(rnorm(prod(nbNodes)), nbNodes[1], nbNodes[2])
# covar2 <- matrix(rnorm(prod(nbNodes)), nbNodes[1], nbNodes[2])
# covarList <- list(covar1 = covar1, covar2 = covar2)

#------------------------------------------------------------------------------------------
#---------------------------- TEST 1  'Bernoulli' model, with nodes covariates
#------------------------------------------------------------------------------------------

test_that("BipartiteSBM_fit 'Bernoulli' model, with nodes covariates (row and col)", {
    ## BIPARTITE BERNOULLI SBM WITH NODES COVARIATES
    means <- matrix(c(
        0.05, 0.95, 0.40,
        0.75, 0.15, 0.60,
        0.30, 0.85, 0.10
    ), 3, 3, byrow = TRUE)
    connectParam <- list(mean = means)

    ## Basic construction with nodes covariates
    mySampler <- BipartiteSBM$new("bernoulli", nbNodes= nbNodes, blockProp = vector('list',2), connectParam=connectParam,
                                  nodesCovarList = nodesCovarList, nodesCovarParam = nodesCovariatesParam
    )
    mySampler$rMemberships(store = TRUE)
    mySampler$rEdges(store = TRUE)

    ## Construction
    mySBM <- BipartiteSBM_fit$new(mySampler$networkData, "bernoulli",
                                  nodesCovarList = nodesCovarList
    )

    ## Checking class
    expect_true(inherits(mySBM, "SBM"))
    expect_true(inherits(mySBM, "BipartiteSBM"))
    expect_true(inherits(mySBM, "BipartiteSBM_fit"))

    ## Checking field access prior to estimation
    expect_equal(mySBM$modelName, "bernoulli")
    expect_equal(unname(mySBM$nbNodes), nbNodes)
    expect_equal(mySBM$dimLabels, c(row = "row", col = "col"))

    ## nodes covariates
    expect_equal(unname(mySBM$nbNodesCovariates[1]), nbNodesCovarRow)
    expect_equal(unname(mySBM$nbNodesCovariates[2]), nbNodesCovarCol)
    expect_equal(dim(mySBM$nodesCovariates[[1]]), c(nbNodes[1], nbNodesCovarRow))
    expect_equal(dim(mySBM$nodesCovariates[[2]]), c(nbNodes[2], nbNodesCovarCol))

    ## Estimation
    BM_out <- mySBM$optimize(estimOptions = list(verbosity = 0, fast = TRUE))
    mySBM$setModel(4)

    ## Expectation
    expect_equal(dim(mySBM$expectation), nbNodes)
    expect_true(all(mySBM$expectation >= 0, na.rm = TRUE))
    expect_true(all(mySBM$expectation <= 1, na.rm = TRUE))

    ## blocks
    expect_equal(names(mySBM$nbBlocks), c("row", "col"))
    expect_true(all(mySBM$nbBlocks >= 1))

    ## nodes covariate parameters
    expect_true(!is.null(mySBM$nodesCovariatesParam))
    expect_true(length(mySBM$nodesCovariatesParam) > 0)
})

#------------------------------------------------------------------------------------------
#---------------------------- TEST 2  'Bernoulli' model, with nodes covariates only on row
#------------------------------------------------------------------------------------------


test_that("BipartiteSBM_fit 'Bernoulli' model, with nodes covariates (row only)", {
    ## BIPARTITE BERNOULLI SBM WITH NODES COVARIATES ON ROW ONLY
    means <- matrix(c(
        0.05, 0.95, 0.40,
        0.75, 0.15, 0.60,
        0.30, 0.85, 0.10
    ), 3, 3, byrow = TRUE)
    connectParam <- list(mean = means)

    nodesCovarRowOnly <- list(row = nodesCovarRow,col=matrix(0,0,0))
    blockPropRowOnly <- vector('list',2)
    blockPropRowOnly[[1]] <- c()
    blockPropRowOnly[[2]] <- blockProp[[2]]
    nodesCovarParam2 <- nodesCovariatesParam
    nodesCovarParam2[[2]] <- c()

    ## Basic construction
    mySampler <- BipartiteSBM$new("bernoulli", nbNodes=nbNodes, blockProp = blockPropRowOnly,  connectParam = connectParam,
                                  nodesCovarList  = nodesCovarRowOnly, nodesCovarParam =nodesCovarParam2
    )
    mySampler$rMemberships(store = TRUE)
    mySampler$rEdges(store = TRUE)

    ## Construction

    mySBM <- BipartiteSBM_fit$new(mySampler$networkData, "bernoulli",
                                  nodesCovarList = nodesCovarRowOnly
    )

    ## Checking class
    expect_true(inherits(mySBM, "BipartiteSBM_fit"))

    ## nodes covariates - only row nodes have covariates
    expect_equal(unname(mySBM$nbNodesCovariates["row"]), nbNodesCovarRow)
    expect_equal(unname(mySBM$nbNodesCovariates["col"]), 0)
    expect_equal(dim(mySBM$nodesCovariates[["row"]]), c(nbNodes[1], nbNodesCovarRow))

    ## Estimation
    BM_out <- mySBM$optimize(estimOptions = list(verbosity = 0, fast = TRUE))
    mySBM$setModel(2)

    ## Expectation
    expect_equal(dim(mySBM$expectation), nbNodes)
})

test_that("BipartiteSBM_fit 'Poisson' model, with nodes covariates", {
    ## BIPARTITE POISSON SBM WITH NODES COVARIATES
    means <- matrix(c(
        10, 5, 7,
        15, 20, 8,
        12, 6, 18
    ), 3, 3, byrow = TRUE)
    connectParam <- list(mean = means)

    ## Basic construction
    mySampler <- BipartiteSBM$new("poisson", nbNodes, blockProp, connectParam,
                                  nodesCovarList = nodesCovarList, nodesCovarParam = nodesCovariatesParam
    )
    mySampler$rMemberships(store = TRUE)
    mySampler$rEdges(store = TRUE)

    ## Construction
    mySBM <- BipartiteSBM_fit$new(mySampler$networkData, "poisson",
                                  nodesCovarList = nodesCovarList
    )

    ## Checking class
    expect_true(inherits(mySBM, "BipartiteSBM_fit"))

    ## nodes covariates
    expect_true(all(mySBM$nbNodesCovariates > 0))
    expect_equal(length(mySBM$nodesCovariates), 2)

    ## Estimation
    BM_out <- mySBM$optimize(estimOptions = list(verbosity = 0, fast = TRUE))
    mySBM$setModel(2)

    ## Expectation
    expect_equal(dim(mySBM$expectation), nbNodes)
    expect_true(all(mySBM$expectation >= 0, na.rm = TRUE))
})

test_that("BipartiteSBM_fit 'Gaussian' model, with nodes covariates", {
    ## BIPARTITE GAUSSIAN SBM WITH NODES COVARIATES
    means <- matrix(c(
        0.5, -0.5, 0.3,
        0.7, -0.2, 0.6,
        0.1, 0.4, -0.3
    ), 3, 3, byrow = TRUE)
    connectParam <- list(mean = means, var = matrix(0.1,3,3))

    ## Basic construction
    mySampler <- BipartiteSBM$new("gaussian", nbNodes, blockProp, connectParam,
                                  nodesCovarList = nodesCovarList, nodesCovarParam=nodesCovariatesParam)
    mySampler$rMemberships(store = TRUE)
    mySampler$rEdges(store = TRUE)

    ## Construction
    mySBM <- BipartiteSBM_fit$new(mySampler$networkData, "gaussian",
                                  nodesCovarList =  nodesCovarList
    )

    ## Checking class
    expect_true(inherits(mySBM, "BipartiteSBM_fit"))

    ## nodes covariates
    expect_equal(unname(mySBM$nbNodesCovariates["row"]), nbNodesCovarRow)
    expect_equal(unname(mySBM$nbNodesCovariates["col"]), nbNodesCovarCol)

    ## Estimation
    BM_out <- mySBM$optimize(estimOptions = list(verbosity = 0, fast = TRUE))
    mySBM$setModel(2)

    ## Expectation
    expect_equal(dim(mySBM$expectation), nbNodes)
})

test_that("BipartiteSBM_fit without nodes covariates returns zero dimensions", {
    ## BIPARTITE BERNOULLI SBM WITHOUT NODES COVARIATES
    means <- matrix(c(
        0.05, 0.95, 0.40,
        0.75, 0.15, 0.60,
        0.30, 0.85, 0.10
    ), 3, 3, byrow = TRUE)
    connectParam <- list(mean = means)

    ## Basic construction without nodes covariates
    mySampler <- BipartiteSBM$new("bernoulli", nbNodes, blockProp, connectParam)
    mySampler$rMemberships(store = TRUE)
    mySampler$rEdges(store = TRUE)

    ## Construction without nodes covariates
    mySBM <- BipartiteSBM_fit$new(mySampler$networkData, "bernoulli")

    ## Check that nbNodesCovariates returns named vector of zeros
    expect_equal(unname(mySBM$nbNodesCovariates["row"]), 0)
    expect_equal(unname(mySBM$nbNodesCovariates["col"]), 0)
    expect_equal(names(mySBM$nbNodesCovariates), c("row", "col"))

    ## nodesCovariates should be empty
    expect_equal(length(mySBM$nodesCovariates), 0)
})
