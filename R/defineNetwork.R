#' Define a network
#'
#' @inheritParams estimateSimpleSBM
#' @param type Type of the matrix, choice between 'simple' and 'bipartite'
#' @param dimLabels an optional list of labels for each dimension (in row, in column)
#' @return a network ready to be part of a multipartite network
#' @examples
#' A <- matrix(rbinom(100,1,.2),10,10)
#' type <- "simple"
#' defineNetwork(A,"poisson",type,directed=TRUE,dimLabels=list("Actor","Actor"))
#' @export
defineNetwork = function(netMat,
                         model        = 'bernoulli',
                         type,
                         directed     = !isSymmetric(netMat),
                         dimLabels    = list(row = "rowLabel", col = "colLabels"),
                         covariates   = list())
  {
  if (!type %in% c("simple","bipartite")) {stop("not allowed type")}

  if (type=="simple")
    mySBM <- SimpleSBM_fit$new(netMat, model, directed, dimLabels, covariates)
   else {mySBM <- SimpleSBM_fit$new(netMat, model, dimLabels, covariates)}

  mySBM
  }

