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

check_boundaries <- function(x, zero = .Machine$double.eps) {
  x[is.nan(x)] <- zero
  x[x > 1 - zero] <- 1 - zero
  x[x <     zero] <-     zero
  x
}

#' @importFrom corrplot corrplot
.prepare_plot <-  function(mat, cl, ordered) {
  if (is.list(cl)) {
    Z1 <- as_indicator(as.factor(cl[[1]]))
    Z2 <- as_indicator(as.factor(cl[[2]]))
    colors <- matrix(-(ncol(Z1) + ncol(Z2)), ncol(Z1), ncol(2));
    colorMat <- Z1 %*% colors %*% t(Z2)
    if (ordered) {
      colorMap <- colorMat[order(cl[[1]]),order(cl[[2]])]
      mat <- mat[order(cl[[1]]), order(cl[[2]])]
    }
  } else {
    Z <- as_indicator(as.factor(cl))
    colors <- matrix(-ncol(Z), ncol(Z), ncol(Z)); diag(colors) <- floor(ncol(Z)/2) + (1:ncol(Z)) # discriminate intra/inter cols
    colorMat <- Z %*% colors %*% t(Z)
    if (ordered) {
      colorMap <- colorMat[order(cl),order(cl)]
      mat <- mat[order(cl), order(cl)]
    }
  }
  mat * colorMap
}

