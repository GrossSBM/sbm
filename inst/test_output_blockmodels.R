## SBM

## --------------------------------------------------------------
## BINARY

## generation of one SBM network
nbNodes  <- 90
blockProp <- c(.5, .25, .25)      # group proportions
connectProb <- diag(.4, 3) + 0.05 # connectivity matrix: affiliation network

## Sampling
mySampler <- sampleSimpleSBM(nbNodes, blockProp, list(mu = connectProb))
adjacencyMatrix <- mySampler$netMatrix

## estimation
my_model <- BM_bernoulli("SBM", adjacencyMatrix, explore_max = 3, explore_max = 3, plotting = "")
my_model$estimate()
my_model_bernoulli <- my_model

npc <- 30 # nodes per class
Q <- 3 # classes
n <- npc * Q # nodes
sigmo <- function(x){1/(1+exp(-x))}
Z<-diag(Q)%x%matrix(1,npc,1)
Mg<-8*matrix(runif(Q*Q),Q,Q)-4
Y1 <- matrix(runif(n*n),n,n)-.5
Y2 <- matrix(runif(n*n),n,n)-.5
M_in_expectation<-sigmo(Z%*%Mg%*%t(Z) + 5*Y1-3*Y2)
M<-1*(matrix(runif(n*n),n,n)<M_in_expectation)
covar_bernoulli <- list(Y1, Y2)
## estimation
my_model <- BM_bernoulli_covariates_fast("SBM", M, list(Y1,Y2), explore_max = 3, explore_max = 3, plotting = "")
my_model$estimate()
my_model_bernoulli_covariates <- my_model

## --------------------------------------------------------------
## POISSON
npc <- 30 # nodes per class
Q <- 3 # classes
n <- npc * Q # nodes
Z<-diag(Q)%x%matrix(1,npc,1)
L<-70*matrix(runif(Q*Q),Q,Q)
M_in_expectation<-Z%*%L%*%t(Z)
M<-matrix(
rpois(
length(as.vector(M_in_expectation)),
as.vector(M_in_expectation))
,n,n)
## estimation
my_model <- BM_poisson("SBM",M,  explore_max = 3, explore_max = 3, plotting = "")
my_model$estimate()
my_model_poisson <- my_model

npc <- 30 # nodes per class
Q <- 3 # classes
n <- npc * Q # nodes
Z<-diag(Q)%x%matrix(1,npc,1)
L<-70*matrix(runif(Q*Q),Q,Q)
M_in_expectation_without_covariates<-Z%*%L%*%t(Z)
Y1 <- matrix(runif(n*n),n,n)
Y2 <- matrix(runif(n*n),n,n)
covar_poisson <- list(Y1, Y2)
M_in_expectation<-M_in_expectation_without_covariates*exp(4.2*Y1-1.2*Y2)
M<-matrix(
rpois(
length(as.vector(M_in_expectation)),
as.vector(M_in_expectation))
,n,n)
## estimation
my_model <- BM_poisson_covariates("SBM",M,list(Y1,Y2),  explore_max = 3, explore_max = 3, plotting = "")
my_model$estimate()
my_model_poisson_covariates <- my_model

## --------------------------------------------------------------
## GAUSSIAN
npc <- 30 # nodes per class
Q <- 3 # classes
n <- npc * Q # nodes
Z<-diag(Q)%x%matrix(1,npc,1)
Mu<-20*matrix(runif(Q*Q),Q,Q)
M<-matrix(rnorm(n*n,sd=10),n,n)+Z%*%Mu%*%t(Z) ## adjacency matrix
## estimation
my_model <- BM_gaussian("SBM",M, explore_max = 3, explore_max = 3, plotting = "" )
my_model$estimate()
my_model_gaussian <- my_model

npc <- 30 # nodes per class
Q <- 3 # classes
n <- npc * Q # nodes
Z<-diag(Q)%x%matrix(1,npc,1)
Mu<-20*matrix(runif(Q*Q),Q,Q)
Y1 <- matrix(runif(n*n),n,n)
Y2 <- matrix(runif(n*n),n,n)
covar_gaussian <- list(Y1, Y2)
M<-matrix(rnorm(n*n,sd=5),n,n)+Z%*%Mu%*%t(Z)+4.2*Y1-1.6*Y2 ## adjacency matrix
## estimation
my_model <- BM_gaussian_covariates("SBM",M,list(Y1,Y2), explore_max = 3, explore_max = 3, plotting = "" )
my_model$estimate()
my_model_gaussian_covariates <- my_model


## BERNOULLI WITH/WITHOUT COVARIATES
theta <- list(pi = my_model_bernoulli$model_parameters[[3]]$pi)
beta <- my_model_bernoulli$model_parameters[[3]]$beta
Z <- my_model_bernoulli$memberships[[3]]$Z
E_bernoulli <- Z %*% theta$pi %*% t(Z)
sum((E_bernoulli - my_model_bernoulli$prediction(3))^2)

beta  <- my_model_bernoulli_covariates$model_parameters[[3]]$beta
theta <- list(m = my_model_bernoulli_covariates$model_parameters[[3]]$m)
Z <- my_model_bernoulli_covariates$memberships[[3]]$Z
E_bernoulli_covariates <- .logistic( Z %*% theta$m %*% t(Z) + GSBM:::roundProduct(simplify2array(covar_bernoulli), beta))
sum((E_bernoulli_covariates - my_model_bernoulli_covariates$prediction(3))^2)

E_bernoulli_covariates <- .logistic(.logit( Z %*% .logistic(theta$m) %*% t(Z)) + GSBM:::roundProduct(simplify2array(covar_bernoulli), beta))

## POISSON WITH/WITHOUT COVARIATES
theta <- list(lambda = my_model_poisson$model_parameters[[3]]$lambda)
beta <- my_model_poisson$model_parameters[[3]]$beta
Z <- my_model_poisson$memberships[[3]]$Z
E_poisson <- Z %*% theta$lambda %*% t(Z)
sum((E_poisson - my_model_poisson$prediction(3))^2)

theta <- list(lambda = my_model_poisson_covariates$model_parameters[[3]]$lambda)
beta <- my_model_poisson_covariates$model_parameters[[3]]$beta
Z <- my_model_poisson_covariates$memberships[[3]]$Z
E_poisson_covariates <- exp( log(Z %*% theta$lambda %*% t(Z))  +  GSBM:::roundProduct(simplify2array(covar_poisson), beta))
sum((E_poisson_covariates - my_model_poisson_covariates$prediction(3))^2)

theta <- list(mu = my_model_gaussian$model_parameters[[3]]$m, sigma2 = my_model_gaussian$model_parameters[[3]]$sigma2)
beta <- my_model_gaussian$model_parameters[[3]]$beta
Z <- my_model_gaussian$memberships[[3]]$Z
E_gaussian <- Z %*% theta$mu %*% t(Z)
sum((E_gaussian - my_model_gaussian$prediction(3))^2)

theta <- list(mu = my_model_gaussian_covariates$model_parameters[[3]]$m, sigma2 = my_model_gaussian_covariates$model_parameters[[3]]$sigma2)
beta <- my_model_gaussian_covariates$model_parameters[[3]]$beta
Z <- my_model_gaussian_covariates$memberships[[3]]$Z
E_gaussian_covariates <- Z %*% theta$mu %*% t(Z) + GSBM:::roundProduct(simplify2array(covar_gaussian), beta)
sum((E_gaussian_covariates - my_model_gaussian_covariates$prediction(3))^2)


