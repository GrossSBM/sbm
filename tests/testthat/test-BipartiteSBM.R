test_that("active binding is working in the class", {
  A <- matrix(rbinom(200,1,.2),20,10)

  myBipartite <- BipartiteSBM$new(incidenceMatrix = A,model = "bernoulli",dimLabels = list("Actor","Stuff"))


  tau1 <- matrix(runif(20*2),20,2)
  tau1 <- tau1 / rowSums(tau1)
  tau2 <- matrix(runif(10*3),10,3)
  tau2 <- tau2 / rowSums(tau2)
  myBipartite$varProb <- list(tau1,tau2)

  myBipartite$blockProp <- list(colMeans(tau1),colMeans(tau2))
  myBipartite$connectParam <- matrix(runif(3*2),3,2)

  expect_equal(myBipartite$dimension,c(20,10))

  expect_equal(myBipartite$memberships[[1]], 1+(tau1[,1]<.5)*1)
  expect_equal(dim(myBipartite$connectParam),c(3,2))
  expect_equal(length(myBipartite$blockProp),2)

})
