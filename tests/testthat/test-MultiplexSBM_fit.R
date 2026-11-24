test_that("Inference for Multiplex networks", {

  if (Sys.info()['sysname'] != "Windows") {

    set.seed(2)
    npc <- 30 # nodes per class
    Q <- 3 # classes
    n <- npc * Q # nodes
    Z<-diag(Q)%x%matrix(1,npc,1)
    P<-matrix(runif(Q*Q),Q,Q)
    A<-1*(matrix(runif(n*n),n,n)<Z%*%P%*%t(Z))
    type <- "simple"
    netA <- defineSBM(A,"bernoulli",type = "simple",directed=TRUE,dimLabels=c("Actor"))
    B <- 1*(matrix(runif(n*n),n,n)<Z%*%P%*%t(Z))
    netB <- defineSBM(B,"bernoulli",type = "simple",dimLabels=c("Actor"))
    myMultiplex <- MultiplexSBM_fit$new(list(netA,netB))
    netC <- defineSBM(B*4,"poisson",type = "simple",dimLabels=c("Actor"))

    expect_equal(myMultiplex$directed, c(TRUE,TRUE))
    expect_equal(myMultiplex$nbNetworks,2)
    expect_equal(myMultiplex$dependentNetwork,FALSE)
    expect_equal(MultiplexSBM_fit$new(list(netA,netB), TRUE)$dependentNetwork,TRUE)
    expect_error(MultiplexSBM_fit$new(list(netA,netC), TRUE))
    expect_equal(MultiplexSBM_fit$new(list(netA,netB,netB), TRUE)$dependentNetwork, TRUE)
    expect_warning(defineSBM(B,"poisson",type = "simple",dimLabels=c("Actor")), "binary")


    currentOptions <- list(
      verbosity     = 0,
      plot          = FALSE,
      nbBlocksRange = list(c(1,10)),
      nbCores       = 2,
      maxiterVE     = 100,
      maxiterVEM    = 100,
      initBM = TRUE
    )

    myMultiplexFitindep <- MultiplexSBM_fit$new(list(netA,netB,netC))
    myMultiplexFitindep$optimize(estimOptions = currentOptions)

    expect_equal(length(myMultiplexFitindep$connectParam),3)

    ## parameters plot: one heat map per layer, the colors of each layer rescaled (different models)
    g <- plot(myMultiplexFitindep, type = "parameters")
    expect_s3_class(g, "ggplot")
    Q <- unname(myMultiplexFitindep$nbBlocks)
    expect_equal(nrow(g$layers[[1]]$data), 3 * Q^2)
    expect_equal(max(g$layers[[1]]$data$fillValue), 1)

    myMultiplexFitdep <- MultiplexSBM_fit$new(list(netA,netB),dependentNet = TRUE)
    currentOptions <- list(
      verbosity     = 0,
      plot          = FALSE,
      explorFactor  = 1.5,
      nbBlocksRange = c(4,Inf),
      nbCores       = 2,
      fast          = TRUE
    )

    myMultiplexFitdep$optimize(estimOptions = currentOptions)
    myMultiplexFitdep$probMemberships

    expect_equal(class(myMultiplexFitdep$memberships),"list")

    expect_equal(length(myMultiplexFitdep$connectParam),4)
    expect_equal(myMultiplexFitdep$dependentNetwork,TRUE)
    expect_s3_class(plot(myMultiplexFitdep, type = "parameters"), "ggplot")


    set.seed(2)
    npc1 <- 30 # nodes per class
    npc2 <- 20
    Q1 <- 2 # classes
    Q2 <- 3
    n1 <- npc1 * Q1 # nodes
    n2 <- npc2 * Q2 # nodes
    Z1 <-diag(Q1)%x%matrix(1,npc1,1)
    Z2 <-diag(Q2)%x%matrix(1,npc2,1)
    P<-matrix(runif(Q1*Q2),Q1,Q2)
    A<-1*(matrix(runif(n1*n2),n1,n2)<Z1%*%P%*%t(Z2))
    netA <- defineSBM(A,"bernoulli",type = "bipartite",directed=TRUE,dimLabels=c("Actor","Object"))
    B <- 1*(matrix(runif(n1*n2),n1,n2)<Z1%*%P%*%t(Z2))
    netB <- defineSBM(B,"bernoulli",type = "bipartite",dimLabels=c("Actor","Object"))
    myMultiplexFitindep <- MultiplexSBM_fit$new(list(netA,netB))
    currentOptions <- list(
      verbosity     = 0,
      nbBlocksRange = list(c(1,10),c(1,10)),
      nbCores       = 2,
      maxiterVE     = 100,
      maxiterVEM    = 100,
      initBM = FALSE
    )
    names(currentOptions$nbBlocksRange) = c("Actor","Object")
    myMultiplexFitindep$optimize(currentOptions)
    expect_equal(length(myMultiplexFitindep$connectParam),2)

  }

})

test_that("Dependent Bernoulli multiplex with more than two layers", {

  if (Sys.info()['sysname'] != "Windows") {

    set.seed(3)
    Q <- 2 ; d <- 3
    combis <- c("000","001","010","011","100","101","110","111")
    ## probabilities of the 8 combinations, for each pair of blocks
    W <- array(rexp(Q * Q * 8), c(Q, Q, 8))
    W[1, 1, ] <- W[1, 1, ] * c(1, 1, 1, 1, 1, 1, 1, 20) # block 1 densely connected in all layers
    W <- sweep(W, 1:2, apply(W, 1:2, sum), "/")
    connectParam <- setNames(lapply(1:8, function(c) W[, , c]), paste0("prob", combis))

    sim <- sampleMultiplexSBM(nbNodes = 60, blockProp = c(.5, .5), nbLayers = d,
                              connectParam = connectParam, model = rep("bernoulli", d),
                              type = "directed", dependent = TRUE)
    expect_length(sim$listSBM, d)

    ## errors on incomplete or non-normalized parameters
    expect_error(sampleMultiplexSBM(nbNodes = 60, blockProp = c(.5, .5), nbLayers = d,
                                    connectParam = connectParam[-1], model = rep("bernoulli", d),
                                    type = "directed", dependent = TRUE), "prob000")
    bad <- connectParam ; bad$prob000 <- bad$prob000 + .1
    expect_error(sampleMultiplexSBM(nbNodes = 60, blockProp = c(.5, .5), nbLayers = d,
                                    connectParam = bad, model = rep("bernoulli", d),
                                    type = "directed", dependent = TRUE), "sum to one")

    fit <- estimateMultiplexSBM(sim$listSBM, dependent = TRUE,
                                estimOptions = list(verbosity = 0, plot = FALSE, nbCores = 1,
                                                    nbBlocksRange = c(3, 3)))
    expect_true(fit$dependentNetwork)
    expect_equal(names(fit$connectParam), paste0("prob", combis))
    Qhat <- unname(fit$nbBlocks)
    expect_equal(Qhat, Q)
    ## the probabilities of the combinations sum to one for each pair of blocks
    expect_equal(Reduce(`+`, fit$connectParam), matrix(1, Qhat, Qhat))
    ## marginal connectivity of layer k: sum over the combinations with a 1 in position k
    pred <- suppressWarnings(fit$predict())
    expect_length(pred, d)
    marg2 <- Reduce(`+`, fit$connectParam[substr(combis, 2, 2) == "1"])
    Z <- fit$probMemberships[[1]]
    expect_equal(pred[[2]], Z %*% marg2 %*% t(Z))
    expect_s3_class(plot(fit, type = "parameters"), "ggplot")
    expect_s3_class(plot(fit, type = "graphon"), "ggplot")
  }
})

test_that("setModel in the dependent multiplex case", {

  if (Sys.info()['sysname'] != "Windows") {
    set.seed(4)
    Q <- 2
    W <- array(rexp(Q * Q * 4), c(Q, Q, 4)); W[1, 1, 4] <- 10 * W[1, 1, 4]
    W <- sweep(W, 1:2, apply(W, 1:2, sum), "/")
    connectParam <- setNames(lapply(1:4, function(c) W[, , c]), c("prob00", "prob01", "prob10", "prob11"))
    sim <- sampleMultiplexSBM(nbNodes = 50, blockProp = c(.5, .5), nbLayers = 2,
                              connectParam = connectParam, model = rep("bernoulli", 2),
                              type = "directed", dependent = TRUE)
    fit <- estimateMultiplexSBM(sim$listSBM, dependent = TRUE,
                                estimOptions = list(verbosity = 0, plot = FALSE, nbCores = 1,
                                                    nbBlocksRange = c(4, 4)))
    best <- unname(fit$nbBlocks)
    fit$setModel(4)
    expect_equal(unname(fit$nbBlocks), 4)
    expect_equal(dim(fit$connectParam$prob11), c(4, 4))
    expect_equal(fit$ICL, fit$storedModels$ICL[fit$storedModels$indexModel == 4])
    expect_length(unique(fit$memberships[[1]]), 4)
    fit$setModel(best)
    expect_equal(unname(fit$nbBlocks), best)
    expect_error(fit$setModel(99), "indexModel")
  }
})
