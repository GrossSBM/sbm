test_that("initializing Multipartite SBM works", {

  set.seed(2)
  #A <- matrix(rbinom(100,1,.2),10,10)
  npc <- 30 # nodes per class
  Q <- 3 # classes
  n <- npc * Q # nodes
  Z<-diag(Q)%x%matrix(1,npc,1)
  P<-matrix(runif(Q*Q),Q,Q)
  A<-1*(matrix(runif(n*n),n,n)<Z%*%P%*%t(Z))
  type <- "simple"
  netA <- defineSBM(A,"bernoulli",type,directed=TRUE,dimLabels=list("Actor","Actor"))
  B <- matrix(rpois(npc*Q*20,2),npc*Q,20)
  type <- "bipartite"
  netB <- defineSBM(B,"poisson",type,directed=TRUE,dimLabels=list("Actor","Stuff"))


  myMBM <- MultipartiteSBM_fit$new(list(netA,netB))
  expect_equal(myMBM$getBM(1)$dimension,Q*c(npc,npc))
  expect_equal(myMBM$getBM(2)$dimension,c(Q*npc,20))

})
