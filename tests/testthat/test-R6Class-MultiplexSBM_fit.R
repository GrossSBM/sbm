test_that("Inference for Multiplex networks", {

  set.seed(2)
  npc <- 30 # nodes per class
  Q <- 3 # classes
  n <- npc * Q # nodes
  Z<-diag(Q)%x%matrix(1,npc,1)
  P<-matrix(runif(Q*Q),Q,Q)
  A<-1*(matrix(runif(n*n),n,n)<Z%*%P%*%t(Z))
  type <- "simple"
  netA <- defineSBM(A,"bernoulli",type = "simple",directed=TRUE,dimLabels=list("Actor","Actor"))
  B <- 1*(matrix(runif(n*n),n,n)<Z%*%P%*%t(Z))
  netB <- defineSBM(B,"bernoulli",type = "simple",dimLabels=list("Actor","Actor"))
  myMultiplex <- MultiplexSBM$new(list(netA,netB))
  netC <- defineSBM(B,"poisson",type = "simple",dimLabels=list("Actor","Actor"))

  myMultiplexFit <- MultiplexSBM_fit$new(list(netA,netB))
  myMultiplexFit$directed

  currentOptions <- list(
    verbosity     = 1,
    nbBlocksRange = lapply(1:myMSBM$nbLabels,function(l){c(1,10)}),
    nbCores       = 2,
    maxiterVE     = 100,
    maxiterVEM    = 100,
    initBM = TRUE
  )

  names(currentOptions$nbBlocksRange) <- myMSBM$dimLabels
  ## Current options are default expect for those passed by the user
  currentOptions[names(estimOptions)] <- estimOptions


  myMultiplexFit$optimize(currentOptions)


  expect_equal(2 * 2, 4)
})
