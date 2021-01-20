test_that("multiplex sampler works", {

  ## Independent multiplex
  Nnodes <- 40
  blockProp <- c(.4,.6)
  nbLayers <- 2
  connectParam <- list(list(mean=matrix(rbeta(4,.5,.5),2,2)),list(mean=matrix(rexp(4,.5),2,2)))
  model <- c("bernoulli","poisson")
  type <- "directed"
  sampMultiplexIndep <- SampleMultiplexSBM(nbNodes = Nnodes,blockProp = blockProp,nbLayers = nbLayers,connectParam = connectParam,model=model,type=type)

  expect_equal(dim(sampMultiplexIndep$listSBM[[2]]$netMatrix),c(Nnodes,Nnodes))
  expect_equal(length(sampMultiplexIndep$memberships[[1]]),Nnodes)


  ## Independent bipartite multiplex
  Nnodes <- c(40,30)
  blockProp <- list(c(.4,.6),rep(.5,2))
  nbLayers <- 2
  connectParam <- list(list(mean=matrix(rbeta(4,.5,.5),2,2)),list(mean=matrix(rexp(4,.5),2,2)))
  model <- c("bernoulli","poisson")
  type <- "bipartite"
  sampMultiplexIndep <- SampleMultiplexSBM(nbNodes = Nnodes,blockProp = blockProp,nbLayers = nbLayers,connectParam = connectParam,model=model,type=type)

  expect_equal(dim(sampMultiplexIndep$listSBM[[2]]$netMatrix),Nnodes)
  expect_equal(length(sampMultiplexIndep$memberships[[1]]),Nnodes[1])
  expect_equal(length(sampMultiplexIndep$memberships[[2]]),Nnodes[2])


  # Dependent bernoulli multiplex
  Q <- 2
  P00<-matrix(runif(Q*Q),Q,Q)
  P10<-matrix(runif(Q*Q),Q,Q)
  P01<-matrix(runif(Q*Q),Q,Q)
  P11<-matrix(runif(Q*Q),Q,Q)
  SumP<-P00+P10+P01+P11
  P00<-P00/SumP
  P01<-P01/SumP
  P10<-P10/SumP
  P11<-P11/SumP
  connectParam$prob00 = P00
  connectParam$prob01 = P01
  connectParam$prob10 = P10
  connectParam$prob11 = P11


})
