#' Define a network
#'
#' @inheritParams estimateSimpleSBM
#' @param type Type of the matrix, choice between 'simple' and 'bipartite'
#' @param dimLabels an optional list of labels for each dimension (in row, in column)
#' @return an object SimpleSBM or BipartiteSBM with the informations required to define a future multipartite network
#' @examples
#' A <- matrix(rbinom(100,1,.2),10,10)
#' type <- "simple"
#' myNet <- defineSBM(A,"poisson",type,directed=TRUE,dimLabels=list("Actor","Actor"))
#' @export
defineSBM = function(netMat,
                     model      = 'bernoulli',
                     type,
                     directed   = !isSymmetric(netMat),
                     dimLabels  = list(row = "rowLabel", col = "colLabel"),
                     covariates = list())
  {
  if (!type %in% c("simple","bipartite")) {stop("type not allowed")}

  #------------ add rownames and colnames if Needed
  if(is.null(rownames(netMat))){rownames(netMat)<- 1:nrow(netMat)}
  if(is.null(colnames(netMat))){colnames(netMat)<- 1:ncol(netMat)}

  if (type == "simple")
    mySBM <- SimpleSBM$new(netMat, model, directed, dimLabels, covariates)
  else
    mySBM <- BipartiteSBM$new(netMat, model, dimLabels, covariates) # SimpleSBM_autre (sampler ou nouveau ?)

  mySBM
  }

