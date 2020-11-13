test_that("optimize for multipartite SBM runs GREMLIN", {

  set.seed(2)
   #A <- matrix(rbinom(100,1,.2),10,10)
  npc <- 30 # nodes per class
  Q <- 3 # classes
  n <- npc * Q # nodes
  Z <- diag(Q)%x%matrix(1,npc,1)
  P <- matrix(runif(Q*Q),Q,Q)
  A <- 1*(matrix(runif(n*n),n,n)<Z%*%P%*%t(Z))
  type <- "simple"
  netA <- defineSBM(A,"bernoulli",type,directed=TRUE,dimLabels=list("Actor","Actor"))
  B <- matrix(rpois(npc*Q*20,2),npc*Q,20)
  type <- "bipartite"
  netB <- defineSBM(B,"poisson",type,directed=TRUE,dimLabels=list("Actor","Stuff"))

  estimOptions = list(initBM = FALSE,verbosity = 0,nbCores = 2)
  Estim  <- estimateMultipartiteSBM(list(netA,netB),estimOptions)



  # private
  #print(E$GREMLINobject)
  expect_equal(length(Estim$getBM(1)$memberships),npc*Q)
  expect_equal(is.list(Estim$getBM(2)$memberships),TRUE)
  expect_equal(length(Estim$getBM(1)$blockProp),length(unique(Estim$getBM(1)$memberships)))
  expect_equal(Estim$getBM(1)$blockProp,Estim$getBM(2)$blockProp[[1]])
  expect_equal(length(Estim$getBM(1)$blockProp),nrow(Estim$getBM(1)$connectParam$mean))
  expect_equal(ncol(Estim$getBM(1)$connectParam$mean),nrow(Estim$getBM(1)$connectParam$mean))
  muAS <- Estim$getBM(2)$connectParam$mean
  if (is.matrix(muAS)){d2 <- ncol(muAS)} else{d2 <- length(muAS)}

  expect_equal(d2,nrow(Estim$getBM(1)$connectParam$mean))
  expect_equal(lengths(Estim$blockProp),Estim$nbBlocks)
  expect_equal(length(Estim$blockProp),Estim$nbLabels)
  expect_equal(length(Estim$connectParam),Estim$nbNetworks)
  expect_equal(lengths(Estim$memberships),Estim$nbNodes)

})
