test_that("multiplex sampler works", {

  ## Independent multiplex
  Nnodes <- 40
  blockProp <- c(.4,.6)
  nbLayers <- 2
  connectParam <- list(list(mean=matrix(rbeta(4,.5,.5),2,2)),list(mean=matrix(rexp(4,.5),2,2)))
  model <- c("bernoulli","poisson")
  type <- "directed"
  sampMultiplexIndep <- sampleMultiplexSBM(nbNodes = Nnodes,blockProp = blockProp,nbLayers = nbLayers,connectParam = connectParam,model=model,type=type)

  expect_equal(dim(sampMultiplexIndep$listSBM[[2]]$networkData),c(Nnodes,Nnodes))
  expect_equal(length(sampMultiplexIndep$memberships[[1]]),Nnodes)


  ## Independent bipartite multiplex
  Nnodes <- c(40,30)
  blockProp <- list(c(.4,.6),rep(.5,2))
  nbLayers <- 2
  connectParam <- list(list(mean=matrix(rbeta(4,.5,.5),2,2)),list(mean=matrix(rexp(4,.5),2,2)))
  model <- c("bernoulli","poisson")
  type <- "bipartite"
  sampMultiplexIndep <- sampleMultiplexSBM(nbNodes = Nnodes,blockProp = blockProp,nbLayers = nbLayers,connectParam = connectParam,model=model,type=type)

  expect_equal(dim(sampMultiplexIndep$listSBM[[2]]$networkData),Nnodes)
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
  connectParam = list()
  connectParam$prob00 = P00
  connectParam$prob01 = P01
  connectParam$prob10 = P10
  connectParam$prob11 = P11
  model = rep("bernoulli",2)
  type = "directed"
  nbLayers = 2
  Nnodes = 40
  blockProp = c(.6,.4)
  sampMultiplexDepBern <- sampleMultiplexSBM(nbNodes = Nnodes,blockProp = blockProp,nbLayers = nbLayers,connectParam = connectParam,model=model,type=type,dependent=TRUE)

  expect_equal(length(sampMultiplexDepBern$memberships[[1]]),Nnodes)
  expect_equal(dim(sampMultiplexDepBern$listSBM[[1]]$networkData),rep(Nnodes,2))

  expect_error(sampleMultiplexSBM(nbNodes = Nnodes,blockProp = blockProp,nbLayers = nbLayers,connectParam = connectParam,model=model,type="undirected",dependent=TRUE))

  Nnodes <- c(40,30)
  blockProp <- list(c(.4,.6),rep(.5,2))
  sampMultiplexDepBern <- sampleMultiplexSBM(nbNodes = Nnodes,blockProp = blockProp,nbLayers = nbLayers,connectParam = connectParam,model=model,type="bipartite",dependent=TRUE)
  expect_equal(length(sampMultiplexDepBern$memberships),2)
  expect_equal(dim(sampMultiplexDepBern$listSBM[[1]]$networkData),Nnodes)


  # dependent Gaussian multiplex
  Q <- 3
  nbLayers <- 2
  connectParam <- list()
  connectParam$mu <- vector("list",nbLayers)
  connectParam$mu[[1]] <- matrix(rnorm(Q*Q),Q,Q)*10
  connectParam$mu[[2]] <- matrix(rnorm(Q*Q),Q,Q)*2
  connectParam$Sigma <- matrix(c(2,1,1,4),nbLayers,nbLayers)
  model <- rep("gaussian",2)
  type <- "directed"
  Nnodes <- 60
  blockProp <- c(.3,.3,.4)
  sampMultiplexDepGau <- sampleMultiplexSBM(nbNodes = Nnodes,blockProp = blockProp,nbLayers = nbLayers,connectParam = connectParam,model=model,type="undirected",dependent=TRUE)

  expect_equal(length(unique(sampMultiplexDepGau$memberships[[1]])),Q)
  expect_equal(sampMultiplexDepGau$listSBM[[1]]$modelName,"gaussian")
  expect_equal(dim(sampMultiplexDepGau$listSBM[[1]]$networkData),rep(Nnodes,2))


  Nnodes <- c(40,30)
  blockProp <- list(c(.4,.6),rep(.5,2))
  Q <- 2
  connectParam$mu[[1]] <- matrix(rnorm(Q*Q),Q,Q)*10
  connectParam$mu[[2]] <- matrix(rnorm(Q*Q),Q,Q)*2
  sampMultiplexDepGau <- sampleMultiplexSBM(nbNodes = Nnodes,blockProp = blockProp,nbLayers = nbLayers,connectParam = connectParam,model=model,type="bipartite",dependent=TRUE)
  expect_equal(length(sampMultiplexDepGau$memberships),2)
  expect_equal(dim(sampMultiplexDepGau$listSBM[[1]]$networkData),Nnodes)
  expect_equal(sampMultiplexDepGau$listSBM[[2]]$modelName,"gaussian")


})
