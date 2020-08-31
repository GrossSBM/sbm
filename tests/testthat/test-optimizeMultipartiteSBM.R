test_that("optimize for multipartite SBM runs GREMLIN", {
  A <- matrix(rbinom(100,1,.2),10,10)
  type <- "simple"
  netA <- defineNetwork(A,"bernoulli",type,directed=TRUE,dimLabels=list("Actor","Actor"))
  B <- matrix(rpois(10*20,2),10,20)
  type <- "bipartite"
  netB <- defineNetwork(B,"poisson",type,directed=TRUE,dimLabels=list("Actor","Stuff"))


  E <- estimateMultipartiteSBM(list(netA,netB))


  expect_equal(2 * 2, 4)
})
