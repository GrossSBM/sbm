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
    if(is.matrix(indicator)){
      cl <- apply(indicator, 1, which.max)
    }
    if(is.vector(indicator)){
      cl <- indicator
    }
  }
  cl
}

##-------------------------------- Some utils function for math
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

#----------------------- RE-ORDERING

order_mbm <- function(list_theta_mean,list_pi,E){

  nbFG <- length(list_pi)
  oAll <-lapply(1:nbFG,FUN = function(f_){
    V <- U <- rep(0,length(list_pi[[f_]]))
    wrow <- which(E[,1]==f_)
    if (length(wrow)>0){U <- c(rowSums(do.call('cbind',lapply(wrow,function(i){list_theta_mean[[i]]%*%list_pi[[E[i,2]]]}))))}
    wcol <- c(which( (E[,2]==f_) & E[,1] != E[,2] ))
    if (length(wcol)>0){V <- c(rowSums(do.call('cbind',lapply(wcol,function(i){c(list_pi[[E[i,1]]]%*% list_theta_mean[[i]])}))))}
    order(U+V,decreasing = TRUE)})
  oAll
}



#-------------------------- Nb of parameters in a MULTIPARTITE model
computeNbConnectParams_MBM <- function(nbBlocks,distrib,E,directed){
  DIR <- directed
  DIR[is.na(DIR)] = TRUE
  nb <- sapply(1:nrow(E),function(i){
    r <- nbBlocks[E[i,1]]*nbBlocks[E[i,2]]
    if (!DIR[i]){r <- r/2 + nbBlocks[E[i,1]]/2}
    if (distrib[i] == 'gaussian'){r <- r*2}
    if (distrib[i] == 'ZIgaussian'){r <- r*3}
    r}
  )
  sum(nb)
}


