
rmse <- function(theta, theta_star) { sqrt(sum((theta - theta_star)^2)/sum(theta_star^2)) }

test_that("initializing Multipartite SBM works", {

  set.seed(2)
  npc <- 30 # nodes per class
  Q <- 3 # classes
  n <- npc * Q # nodes
  Z <- diag(Q)%x%matrix(1,npc,1)
  P <- matrix(runif(Q*Q),Q,Q)
  A <- 1*(matrix(runif(n*n),n,n)<Z%*%P%*%t(Z))
  B <- matrix(rpois(npc*Q*20,2),npc*Q,20)

  netA <- defineSBM(A, "bernoulli", type = "simple"    , directed=TRUE, dimLabels = list(row = "Actor", col = "Actor"))
  netB <- defineSBM(B, "poisson"  , type = "bipartite", dimLabels = list(row = "Actor", col = "Stuff"))
  myMBM <- MultipartiteSBM_fit$new(list(netA,netB))

  ## Checking class
  expect_true(inherits(myMBM, "MultipartiteSBM"))
  expect_true(inherits(myMBM, "MultipartiteSBM_fit"))

  ## Checking field access and format prior to estimation
  ## parameters
  expect_equal(myMBM$modelName, c('bernoulli', 'poisson'))
  expect_true(is.character(myMBM$modelName))
  expect_equal(unname(myMBM$dimension), c(Q*npc,20))
  expect_equal(unname(myMBM$nbNodes) , c(Q*npc,20))
  expect_equal(myMBM$directed, c(TRUE,NA))
  expect_equal(myMBM$nbNetworks,2)
  expect_equal(myMBM$networkList[[1]]$dimension,Q*c(npc,npc))
  expect_equal(myMBM$networkList[[2]]$dimension,c(Q*npc,20))
  expect_equal(unname(myMBM$architecture), matrix(c(1,1,1,2), 2,2))
  expect_equivalent(myMBM$blockProp, list(NULL, NULL))
  expect_equivalent(myMBM$connectParam, list(list(mean = matrix(NA)), list(mean = matrix(NA))))

  # S3 methods
  expect_silent(plot(myMBM, type = "data"))
  expect_equal(coef(myMBM, 'connectivity'), myMBM$connectParam)
  expect_equal(coef(myMBM, 'block')       , myMBM$blockProp)

  ## Estimation-----------------------------------------------------------------
  estimOptions = list(initBM = FALSE,verbosity = 0,nbCores = 2)
  myMBM$optimize(estimOptions)

  ## Field set after optimization
  expect_equal(length(myMBM$networkList[[1]]$memberships), npc*Q)
  expect_equal(is.list(myMBM$networkList[[2]]$memberships),TRUE)
  expect_equal(length(myMBM$networkList[[1]]$blockProp),
               length(unique(myMBM$networkList[[1]]$memberships)))
  expect_equal(myMBM$networkList[[1]]$blockProp,
               myMBM$networkList[[2]]$blockProp[[1]])
  expect_equal(length(myMBM$networkList[[1]]$blockProp),
               nrow(myMBM$networkList[[1]]$connectParam$mean))
  expect_equal(ncol(myMBM$networkList[[1]]$connectParam$mean),
               nrow(myMBM$networkList[[1]]$connectParam$mean))

  muAS <- myMBM$networkList[[2]]$connectParam$mean
  expect_equal(ifelse(is.matrix(muAS), nrow(muAS), length(muAS)),
               nrow(myMBM$networkList[[1]]$connectParam$mean))
  expect_equal(lengths(myMBM$blockProp), myMBM$nbBlocks)
  expect_equal(length(myMBM$blockProp), myMBM$nbLabels)
  expect_equal(length(myMBM$connectParam), myMBM$nbNetworks)
  expect_equal(lengths(myMBM$memberships), myMBM$nbNodes)
  expect_lt(myMBM$loglik, 0)
  expect_lt(myMBM$ICL, 0)
  expect_lt(myMBM$ICL, myMBM$loglik)

  # S3 methods
  expect_silent(plot(myMBM, type = "data"))
  expect_silent(plot(myMBM, type = "meso"))
  expect_silent(plot(myMBM, type = "expected"))
  expect_equal(coef(myMBM, 'connectivity'), myMBM$connectParam)
  expect_equal(coef(myMBM, 'block')       , myMBM$blockProp)

  ## correctness
  expect_lt(rmse(myMBM$connectParam[[1]]$mean, netA$connectParam$mean), 0.01)
  expect_lt(rmse(myMBM$connectParam[[2]]$mean, netB$connectParam$mean), 0.01)
  expect_lt(1 - aricode::ARI(myMBM$memberships$Actor, netA$memberships), 0.05)

})