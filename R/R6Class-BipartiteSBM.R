#' R6 virtual class for SBM representation
#'
#' @import R6
#' @include R6Class-SBM.R
#' @export
BipartiteSBM <- # this class inherit from SBM and allow to use multipartite as a list of SimpleSBM and BipartiteSBM
  R6::R6Class(classname = "BipartiteSBM",
              inherit = SBM,
              private = list(
                Y = NULL,
                tau = NULL

              ),
              public = list(
                #' @description constructor for a Simple SBM fit
                #' @param incidenceMatrix incidence (weighted) matrix
                #' @param model character (\code{'bernoulli'}, \code{'poisson'}, \code{'gaussian'})
                #' @param directed logical, directed network or not. In not, \code{incidenceMatrix} must be symmetric.
                #' @param dimLabels list of labels of each dimension (in row, in columns)
                #' @param covarList and optional list of covariates, each of whom must have the same dimension as \code{incidenceMatrix}
                initialize = function(incidenceMatrix, model,  dimLabels=list(row="rowLabel", col="colLabel"), covarList=list()) {

                  ## INITIALIZE THE SBM OBJECT ACCORDING TO THE DATA
                  super$initialize(model = model, dimension = dim(incidenceMatrix), dimLabels = dimLabels, covarList = covarList)
                  private$Y <- incidenceMatrix

                }
              ),
              active = list(
                #' @field varProb   variational probabilities of nodes being in a block. Either return the probabilities of clustering or may be used to set that probabilites
                  varProb    = function(value) {if (missing(value)) return(private$tau) else {
                  stopifnot(is.list(value))
                  stopifnot(nrow(value[[1]])==private$dim[1])
                  stopifnot(nrow(value[[2]])==private$dim[2])
                  private$tau <- value}},
                #' @field memberships list of size 2: vector of memberships in row, in column.
                memberships = function(value) {lapply(private$tau, as_clustering)}
              )
  )
