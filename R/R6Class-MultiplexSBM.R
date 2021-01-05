#' R6 Class definition of a Multiplex SBM
#'
#' R6 virtual class for Multiplex SBM representation
#'
#' @export
MultiplexSBM <-
  R6::R6Class(classname = "MultiplexSBM",
              inherit = MultipartiteSBM,
              private = list(
                depStructure = NULL
              ),
              public = list(
                #' @description constructor for Multipartite SBM
                #' @param listSBM list of SimpleSBM or BipartiteSBM
                #' @param memberships list of memberships for each node in each function group.Default value is NULL
                #' @param dep boolean indicating whether dependence is assumed between networks beyond the common dependence on the latent variables
                initialize = function(listSBM, memberships=NULL,dep=FALSE) {
                  # SANITY CHECK
                  # check whether the multipartite at hand is actually a multiplex
                  labs <- sapply(listSBM,function(net) net$dimLabels)
                  if (length(unique(labs[1,]))>1 | length(unique(labs[2,]))>1)
                  {
                    print("list of networks provided does not correspond to a Multiplex architecture")
                    stop()
                  }
                  super$initialize(listSBM, memberships)

                   # SANITY CHECK dependence structure
                  if (dep)
                  {
                    if (length(unique(self$directed))>1)
                    {
                      print("in the dependent case, all networks should be either directed or not directed")
                      stop()
                    }

                    dBern = isTRUE(all.equal(self$modelName,rep("bernoulli",self$nbNetworks)))
                    dGauss = isTRUE(all.equal(self$modelName,rep("gaussian",self$nbNetworks)))

                    if (!(dGauss | (dBern&self$nbNetworks==2)))
                    {
                      print("dependency in multiplex network is only handled for Gaussian distribution or a bivariate Bernoulli distribution")
                      stop()
                    }
                  }
                  private$depStructure <- dep

                }
              ),
              active = list(
                modelDependence = function(value) {private$depStructure}
              )
  )



