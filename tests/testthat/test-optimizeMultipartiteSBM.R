test_that("optimize for multipartite SBM runs GREMLIN", {
  A <- matrix(rbinom(100,1,.2),10,10)
  type <- "simple"
  netA <- defineSBM(A,"bernoulli",type,directed=TRUE,dimLabels=list("Actor","Actor"))
  B <- matrix(rpois(10*20,2),10,20)
  type <- "bipartite"
  netB <- defineSBM(B,"poisson",type,directed=TRUE,dimLabels=list("Actor","Stuff"))



  E <- estimateMultipartiteSBM(list(netA,netB))

  # private
  #print(E$GREMLINobject)

  E$getBM(i = 1)$memberships

  expect_equal(2 * 2, 4)
})
