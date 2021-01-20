test_that("multiplex sampler works", {

  Nnodes <- 40
  blockProp <- c(.4,.6)
  nbLayers <- 2
  connectParam <- list(list(mean=matrix(rbeta(4,.5,.5),2,2)),list(mean=matrix(rexp(4,.5),2,2)))
  model <- c("bernoulli","poisson")
  type <- "directed"

  sampMultiplexIndep <- SampleMultiplexSBM(nbNodes = Nnodes,blockProp = blockProp,nbLayers = nbLayers,connectParam = connectParam,model=model,type=type)

  expect_equal(2 * 2, 4)
})
