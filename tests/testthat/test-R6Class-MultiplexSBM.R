test_that("initializing Multipartite SBM works", {

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

  expect_equal(myMultiplex$directed, c(TRUE,TRUE))
  expect_equal(myMultiplex$nbNetworks,2)
  expect_equal(myMultiplex$modelDependence,FALSE)
  expect_equal(MultiplexSBM$new(list(netA,netB),dep=TRUE)$modelDependence,TRUE)
  expect_error(MultiplexSBM$new(list(netA,netC),dep=TRUE))
  expect_error(MultiplexSBM$new(list(netA,netB,netB),dep=TRUE))


})
