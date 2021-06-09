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
    netC <- defineSBM(B,"poisson",type = "simple",dimLabels=c("Actor"))

    expect_equal(myMultiplex$directed, c(TRUE,TRUE))
    expect_equal(myMultiplex$nbNetworks,2)
    expect_equal(myMultiplex$dependentNetwork,FALSE)
    expect_equal(MultiplexSBM_fit$new(list(netA,netB), TRUE)$dependentNetwork,TRUE)
    expect_error(MultiplexSBM_fit$new(list(netA,netC), TRUE))
    expect_error(MultiplexSBM_fit$new(list(netA,netB,netB), TRUE))

    currentOptions <- list(
      verbosity     = 1,
      nbBlocksRange = list(c(1,10)),
      nbCores       = 2,
      maxiterVE     = 100,
      maxiterVEM    = 100,
      initBM = TRUE
    )

    myMultiplexFitindep <- MultiplexSBM_fit$new(list(netA,netB,netC))
    myMultiplexFitindep$optimize(estimOptions = currentOptions)

    expect_equal(length(myMultiplexFitindep$connectParam),3)

    myMultiplexFitdep <- MultiplexSBM_fit$new(list(netA,netB),dependentNet = TRUE)
    currentOptions <- list(
      verbosity     = 3,
      plot          = TRUE,
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
      verbosity     = 1,
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
