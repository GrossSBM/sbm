test_that("Inference for Multiplex networks", {

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

  myMultiplexFitdep <- MultiplexSBM_fit$new(list(netA,netB),dep = TRUE)
  currentOptions <- list(
    verbosity     = 3,
    plot          = TRUE,
    explorFactor  = 1.5,
    nbBlocksRange = c(4,Inf),
    nbCores       = 2,
    fast          = TRUE
  )

  myMultiplexFitdep$optimize(estimOptions = currentOptions)




 expect_equal(length(myMultiplexFitdep$connectParam),4)

 expect_equal(myMultiplexFitdep$dependentNetwork,TRUE)
})
