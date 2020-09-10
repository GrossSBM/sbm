test_that("optimize for multipartite SBM runs GREMLIN", {

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



  E <- estimateMultipartiteSBM(list(netA,netB))

  # private
  #print(E$GREMLINobject)

  expect_equal(length(E$getBM(1)$memberships),npc*Q)
  expect_equal(is.list(E$getBM(2)$memberships),TRUE)

  expect_equal(length(E$getBM(1)$blockProp),length(unique(E$getBM(1)$memberships)))
  expect_equal(E$getBM(1)$blockProp,E$getBM(2)$blockProp[[1]])

  expect_equal(length(E$getBM(1)$blockProp),nrow(E$getBM(1)$connectParam))
  expect_equal(ncol(E$getBM(1)$connectParam),nrow(E$getBM(1)$connectParam))
  expect_equal(nrow(E$getBM(2)$connectParam),nrow(E$getBM(1)$connectParam))

})
