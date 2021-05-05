#' Define a network
#'
#' @inheritParams estimateSimpleSBM
#' @param type Type of the matrix, choice between 'simple' and 'bipartite'
#' @param dimLabels an optional vector of labels for each dimension (in row, in column). Default value = c('row' = row,'col'= col)
#' @return an object SimpleSBM or BipartiteSBM with the informations required to define a future multipartite network
#' @examples
#' A <- matrix(rbinom(100,1,.2), 10, 10)
#' myNet <- defineSBM(A, "poisson", "simple", TRUE, "Actor")
#' @export
defineSBM = function(netMat,
                     model      = 'bernoulli',
                     type       = ifelse(ncol(netMat) == nrow(netMat), "simple", "bipartite"),
                     directed   = !isSymmetric(netMat),
                     dimLabels  = c(row = "row", col = "col"),
                     covariates = list())
  {

  ## sanity check
  if (!type %in% c("simple","bipartite")) {stop("type not allowed")}

  #------------ add rownames and colnames if needed
  if(is.null(rownames(netMat))) rownames(netMat) <- 1:nrow(netMat)
  if(is.null(colnames(netMat))) colnames(netMat) <- 1:ncol(netMat)

  if (type == "simple")
    mySBM <- SimpleSBM_fit$new(netMat, model = model, directed = directed, dimLabels = dimLabels, covarList  = covariates)
  else
    mySBM <- BipartiteSBM_fit$new(netMat, model = model, dimLabels = dimLabels, covarList = covariates)

  mySBM
  }

