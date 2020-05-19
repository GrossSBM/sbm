as_indicator <- function(clustering) {
  K <- length(unique(clustering))
  N  <- length(clustering)
  Z <- matrix(0, N, K)
  Z[cbind(seq.int(N), clustering)] <- 1
  Z
}

as_clustering <- function(indicator) {
  cl <- apply(indicator, 1, which.max)
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


