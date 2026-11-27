## Additional initializations of the exploration of blockmodels (option 'classifInit' of
## estimateSimpleSBM() and estimateBipartiteSBM()).
##
## The default exploration of blockmodels starts from one block and splits / merges blocks, with
## spectral clustering initializations. It can miss partitions that require several blocks to be
## split at the same time (e.g. two compartments in a food web: consumers and their resources must be
## separated together). These functions provide other initial clusterings (k-means on the rows /
## columns of the network, or clusterings supplied by the user), run the VEM algorithm from them,
## and then let blockmodels split and merge again from the new solutions.

## Profiles of the nodes used by k-means: missing values replaced by the mean of the network
.profiles_for_kmeans <- function(X) {
  X[is.na(X)] <- mean(X, na.rm = TRUE)
  X
}

## k-means clustering of the rows of X in K groups (NULL if not possible)
#' @importFrom stats kmeans
.kmeans_classif <- function(X, K, nstart = 20) {
  if (K == 1) return(rep(1L, nrow(X)))
  if (K >= nrow(unique(X))) return(NULL)
  tryCatch(kmeans(X, K, nstart = nstart, iter.max = 50)$cluster, error = function(e) NULL)
}

## Check a clustering supplied by the user
.check_classif <- function(cl, n, what) {
  if (!is.atomic(cl) || length(cl) != n || anyNA(cl))
    stop("classifInit: each initial clustering of the ", what, " must be a vector of length ", n,
         " without missing values.", call. = FALSE)
  as.integer(factor(cl))
}

## Run the VEM algorithm of a blockmodels object from additional initializations, then explore again
## (splits and merges) from the new solutions.
## inits: list of blockmodels membership objects; Q: index of each one in blockmodels (number of
## blocks for SBM, total number of blocks for LBM).
.add_inits_to_BM <- function(BMobject, inits, Q) {
  if (length(inits) == 0) return(invisible(BMobject))
  ## blockmodels stores the models by index Q = 1, 2, ...: make sure all the indexes up to max(Q)
  ## have been explored before adding the new solutions
  if (max(Q) > length(BMobject$ICL)) {
    BMobject$explore_min <- max(Q)
    BMobject$estimate()
  }
  for (q in sort(unique(Q))) BMobject$do_with_inits(inits[Q == q], q, 1)
  BMobject$estimate()
  invisible(BMobject)
}

## Initializations for a simple SBM
.classifInit_simple <- function(classifInit, BMobject, Y, directed) {
  membership <- ifelse(directed, "SBM", "SBM_sym")
  n <- nrow(Y)
  if (identical(classifInit, "kmeans")) {
    Yp <- Y; diag(Yp) <- NA
    X <- .profiles_for_kmeans(if (directed) cbind(Yp, t(Yp)) else Yp)
    Ks <- seq(2, max(2, length(BMobject$ICL)))
    classifs <- lapply(Ks, function(K) .kmeans_classif(X, K))
  } else {
    if (!is.list(classifInit)) classifInit <- list(classifInit)
    classifs <- lapply(classifInit, .check_classif, n = n, what = "nodes")
  }
  classifs <- Filter(Negate(is.null), classifs)
  nodes_covar <- BMobject$nodes_covariates[["node"]]
  inits <- lapply(classifs, function(cl)
    getRefClass(membership)(classif = cl, nodes_covar = nodes_covar))
  .add_inits_to_BM(BMobject, inits, Q = sapply(classifs, max))
}

## Initializations for a bipartite SBM (LBM: blockmodels indexes the models by K_row + K_col)
.classifInit_bipartite <- function(classifInit, BMobject, Y) {
  dims <- dim(Y)
  if (identical(classifInit, "kmeans")) {
    Qmax <- max(3, length(BMobject$ICL))
    Xrow <- .profiles_for_kmeans(Y)
    Xcol <- .profiles_for_kmeans(t(Y))
    rows <- lapply(seq_len(Qmax - 1), function(K) .kmeans_classif(Xrow, K))
    cols <- lapply(seq_len(Qmax - 1), function(K) .kmeans_classif(Xcol, K))
    pairs <- expand.grid(K1 = seq_len(Qmax - 1), K2 = seq_len(Qmax - 1))
    pairs <- pairs[pairs$K1 + pairs$K2 >= 3 & pairs$K1 + pairs$K2 <= Qmax, ]
    classifs <- lapply(seq_len(nrow(pairs)), function(p) {
      r <- rows[[pairs$K1[p]]]; c <- cols[[pairs$K2[p]]]
      if (is.null(r) || is.null(c)) NULL else list(r, c)
    })
  } else {
    if (!is.list(classifInit) || length(classifInit) == 0)
      stop("classifInit must be \"kmeans\", a list (row, col) of two clusterings, or a list of such lists.",
           call. = FALSE)
    if (!is.list(classifInit[[1]])) classifInit <- list(classifInit)
    classifs <- lapply(classifInit, function(cl) {
      if (length(cl) != 2) stop("classifInit: each initial clustering must be a list of two vectors (row, col).",
                                call. = FALSE)
      list(.check_classif(cl[[1]], dims[1], "rows"), .check_classif(cl[[2]], dims[2], "columns"))
    })
  }
  classifs <- Filter(Negate(is.null), classifs)
  inits <- lapply(classifs, function(cl)
    getRefClass("LBM")(classif = cl, row_covar = BMobject$nodes_covariates[["row"]],
                       col_covar = BMobject$nodes_covariates[["col"]]))
  .add_inits_to_BM(BMobject, inits, Q = sapply(classifs, function(cl) max(cl[[1]]) + max(cl[[2]])))
}
