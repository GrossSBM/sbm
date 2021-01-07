test_that("active binding is working in the class", {

  A <- matrix(rbinom(100,1,.2),10,10)

  mySimple <- SimpleSBM$new(adjacencyMatrix = A,model = "bernoulli",directed = TRUE,dimLabels = list("Actor","Actor"))
  p <- runif(10)
  mySimple$varProb <- matrix(c(p,1-p),10,2,byrow=FALSE)
  mySimple$blockProp <- c(.3,.7)
  mySimple$connectParam <- matrix(runif(4),2,2)


  expect_equal(mySimple$memberships, 1+(p<.5)*1)
  expect_equal(dim(mySimple$connectParam),c(2,2))
  expect_equal(length(mySimple$blockProp),2)

})
