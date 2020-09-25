as_indicator <- function(clustering) {
  K <- length(unique(clustering))
  N  <- length(clustering)
  Z <- matrix(0, N, K)
  Z[cbind(seq.int(N), clustering)] <- 1
  Z
}

as_clustering <- function(indicator) {
  if (is.null(indicator)) {
    cl <- numeric(0)
  } else {
    cl <- apply(indicator, 1, which.max)
  }
  cl
}

## Some utils function for math
.xlogx <- function(x) ifelse(x < .Machine$double.eps, 0, x*log(x))
.logistic <- function(x) {1/(1 + exp(-x))}
.logit    <- function(x) {log(x/(1 - x))}

.softmax <- function(x) {
  b <- max(x)
  exp(x - b) / sum(exp(x - b))
}

.na2zero <- function(x) {
  x[is.na(x)] <- 0
  x
}

check_boundaries <- function(x, zero = .Machine$double.eps) {
  x[is.nan(x)] <- zero
  x[x > 1 - zero] <- 1 - zero
  x[x <     zero] <-     zero
  x
}

#----------------------- PREDICTION
predict_sbm <- function(nbNodes,nbCovariates,link,invlink,tau,theta_mean,covarEffect,covarList){


  stopifnot(is.list(covarList), nbCovariates == length(covarList))
  mu <- tau %*% theta_mean %*% t(tau)
  if (nbCovariates > 0) {
    all(sapply(covarList, nrow) == nbNodes, sapply(covarList, ncol) == nbNodes)
    mu <- invlink(link(mu) + covarEffect)
  }
  mu
}


predict_lbm <- function(dimension,nbCovariates,link,invlink,tau,theta_mean,covarEffect,covarList){

    stopifnot(!is.null(tau[[1]]), !is.null(tau[[2]]), !is.null(theta_mean))
    stopifnot(is.list(covarList),  nbCovariates == length(covarList))

    if (length(covarList) > 0) {
      stopifnot(all(sapply(covarList, nrow) == dimension[1]),
                all(sapply(covarList, ncol) == dimension[2]))
    }
    mu <- tau[[1]] %*% theta_mean %*% t(tau[[2]])
    if (length(covarList) > 0) mu <- invlink(link(mu) + covarEffect)
    mu
}

order_sbm <- function(theta_mean,pi){
  o <- order(theta_mean %*% pi, decreasing = TRUE)
  return(o)
}

order_lbm <- function(theta_mean,pi){
  oRow <- order(theta_mean %*% pi[[2]], decreasing = TRUE)
  oCol <- order(pi[[1]] %*% theta_mean, decreasing = TRUE)
  return(list(row  = oRow, col = oCol))
}



#-------------------------- Nb of parameters in a MULTIPARTITE model
computeNbConnectParams_MBM <- function(nbBlocks,distrib,E,directed){
  DIR <- directed
  DIR[is.na(DIR)] = TRUE
  nb <- sapply(1:nrow(E),function(i){
    r <- nbBlocks[E[i,1]]*nbBlocks[E[i,2]]
    if (!DIR[i]){r <- r/2 + nbBlocks[E[i,1]]/2}
    if (distrib[i] == 'gaussian'){r <- r*2}
    r}
  )
  sum(nb)
}


