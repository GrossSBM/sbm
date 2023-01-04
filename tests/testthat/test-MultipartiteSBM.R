
rmse <- function(theta, theta_star) { sqrt(sum((theta - theta_star)^2)/sum(theta_star^2)) }

test_that("initializing Multipartite SBM works", {

  if (Sys.info()['sysname'] != "Windows") {
    set.seed(2)
    npc <- 30 # nodes per class
    Q <- 3 # classes
    n <- npc * Q # nodes
    Z <- diag(Q)%x%matrix(1,npc,1)
    P <- matrix(runif(Q*Q),Q,Q)
    A <- 1*(matrix(runif(n*n),n,n)<Z%*%P%*%t(Z))
    B <- matrix(rpois(npc*Q*20,2),npc*Q,20)

    netA <- defineSBM(A, "bernoulli", type = "simple"   , directed=TRUE, dimLabels = c("Actor"))
    netB <- defineSBM(B, "poisson"  , type = "bipartite", dimLabels = c("Actor", "Stuff"))
    myMBM <- MultipartiteSBM_fit$new(list(netA,netB))

    ## Checking class
    expect_true(inherits(myMBM, "SBM"))
    expect_true(inherits(myMBM, "MultipartiteSBM"))
    expect_true(inherits(myMBM, "MultipartiteSBM_fit"))

    ## Checking field access and format prior to estimation
    ## parameters
    expect_equal(myMBM$modelName, c('bernoulli', 'poisson'))
    expect_true(is.character(myMBM$modelName))
    expect_equal(unname(myMBM$nbNodes) , c(Q*npc,20))
    expect_equal(myMBM$directed, c(TRUE,NA))
    expect_equal(myMBM$nbNetworks,2)
    expect_equal(unname(myMBM$networkData[[1]]$nbNodes),Q*npc)
    expect_equal(unname(myMBM$networkData[[2]]$nbNodes),c(Q*npc,20))
    expect_equal(unname(myMBM$architecture), matrix(c(1,1,1,2), 2,2))
    if (packageVersion("purrr") >= "1.0.0") {
      expect_equal(myMBM$blockProp, list(numeric(0), list(numeric(0), numeric(0))))
    }
    expect_equal(myMBM$connectParam,
                      list(list(mean = matrix(0,0,0)), list(mean = matrix(0,0,0))))
    # S3 methods
    ##  expect_silent(plot(myMBM, type = "data"))
    expect_equal(coef(myMBM, 'connectivity'), myMBM$connectParam)
    expect_equal(coef(myMBM, 'block')       , myMBM$blockProp)

    ## Estimation-----------------------------------------------------------------
    estimOptions = list(initBM = FALSE,verbosity = 0,nbCores = 2)
    myMBM$optimize(estimOptions)

    ## Field set after optimization
    expect_equal(length(myMBM$networkData[[1]]$memberships), npc*Q)
    expect_equal(is.list(myMBM$networkData[[2]]$memberships),TRUE)
    expect_equal(length(myMBM$networkData[[1]]$blockProp),
                 length(unique(myMBM$networkData[[1]]$memberships)))
    expect_equal(myMBM$networkData[[1]]$blockProp,
                 myMBM$networkData[[2]]$blockProp[[1]])
    expect_equal(length(myMBM$networkData[[1]]$blockProp),
                 nrow(myMBM$networkData[[1]]$connectParam$mean))
    expect_equal(ncol(myMBM$networkData[[1]]$connectParam$mean),
                 nrow(myMBM$networkData[[1]]$connectParam$mean))

    muAS <- myMBM$networkData[[2]]$connectParam$mean
    expect_equal(ifelse(is.matrix(muAS), nrow(muAS), length(muAS)),
                 nrow(myMBM$networkData[[1]]$connectParam$mean))
    expect_equal(lengths(myMBM$blockProp), myMBM$nbBlocks)
    expect_equal(length(myMBM$blockProp), length(myMBM$dimLabels))
    expect_equal(length(myMBM$connectParam), myMBM$nbNetworks)
    expect_equal(lengths(myMBM$memberships), myMBM$nbNodes)
    expect_lt(myMBM$loglik, 0)
    expect_lt(myMBM$ICL, 0)
    expect_lt(myMBM$ICL, myMBM$loglik)

    # S3 methods
    ## expect_silent(plot(myMBM, type = "data"))
    expect_silent(plot(myMBM, type = "meso"))
    ##  expect_silent(plot(myMBM, type = "expected"))
    expect_equal(coef(myMBM, 'connectivity'), myMBM$connectParam)
    expect_equal(coef(myMBM, 'block')       , myMBM$blockProp)

    ## correctness
    expect_lt(rmse(myMBM$connectParam[[1]]$mean, netA$connectParam$mean), 0.01)
    expect_lt(rmse(myMBM$connectParam[[2]]$mean, netB$connectParam$mean), 0.01)
    expect_lt(1 - aricode::ARI(myMBM$memberships$Actor, netA$memberships), 0.05)

  }

})
