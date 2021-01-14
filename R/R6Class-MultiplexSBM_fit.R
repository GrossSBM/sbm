#' R6 Class definition of a Multiplex SBM fit
#'
#' This class is designed to give a representation and adjust a Multiplex SBM fitted with GREMLIN.
#'
#' @import R6 GREMLINS
#' @export
MultiplexSBM_fit <-
  R6::R6Class(
    classname = "MultiplexSBM_fit",
    inherit = MultipartiteSBM_fit,
    # fields for internal use (referring to the mathematical notation)
    private = list(
      BMobject = NULL,
      dependent =NULL,
      import_from_BM  = function(index = which.max(private$BMobject$ICL)) {
        private$J     <- private$BMobject$PL[index]
        private$vICL  <- private$BMobject$ICL[index]
        parameters    <- private$BMobject$model_parameters[[index]]
        private$theta <- switch(private$BMobject$model_name,
                "gaussian_multivariate" = list(mean=parameters$mu,cov=parameters$Sigma),
                "bernoulli_multiplex"   = list(prob00=parameters$pi$`00`,prob01=parameters$pi$`01`,
                                               prob00=parameters$pi$`10`,prob10=parameters$pi$`11`)
        )},
      import_from_BM_Simple = function(index = which.max(private$BMobject$ICL)) { # a function updating the Class
        private$import_from_BM(index)
        private$tau <- private$BMobject$memberships[[index]]$Z
        private$pi  <- colMeans(private$tau)
      },
      import_from_BM_Bipartite  = function(index = which.max(private$BMobject$ICL)) {
        private$import_from_BM(index)
        private$tau <- list(
          row = private$BMobject$memberships[[index]]$Z1,
          col = private$BMobject$memberships[[index]]$Z2
        )
        private$pi  <- lapply(private$tau, colMeans)
      }
    ),
    #-----------------------------------------------
    public = list(
      #' @description constructor for Multiplex SBM
      #' @param netList list of SBM object with
      #' @param dependentNet boolean indicating whether dependence is assumed between networks beyond the common dependence on the latent variables
      initialize = function(netList, dependentNet = FALSE) {

        # check whether the multipartite at hand is actually a multiplex
        lab_per_row <- map(netList, "dimLabels") %>% map(~unique(unlist(.x))) %>% map_int(length)
        if (any(lab_per_row > 1))
          stop("list of networks provided does not correspond to a Multiplex architecture")
        super$initialize(netList)

        # CHECKING dependence structure
        if (dependentNet) {
          if (! ( all(self$directed == TRUE) | all(self$directed == FALSE)) )
            stop("in the dependent case, all networks should be either directed or not directed")

          dBern  <- isTRUE(all.equal(self$modelName, rep("bernoulli", self$nbNetworks)))
          dGauss <- isTRUE(all.equal(self$modelName, rep("gaussian" , self$nbNetworks)))
          if (!(dGauss | (dBern&self$nbNetworks == 2)))
            stop("dependency in multiplex network is only handled for Gaussian distribution or a bivariate Bernoulli distribution")
        }
        private$dependent <- dependentNet
      },
      #' @description estimation of multipartiteSBM via GREMLINS
      #' @param estimOptions options for MultipartiteBM
      #' @details The list of parameters \code{estimOptions} essentially tunes the optimization process and the variational EM algorithm, with the following parameters
      #'  \itemize{
      #'  \item{"nbCores"}{integer for number of cores used.  Default is 2}
      #'  \item{"verbosity"}{integer for verbosity (0, 1). Default is 1}
      #'  \item{"nbBlocksRange"}{List of length the number of functional groups, each element supplying the minimal and maximal number of blocks to be explored. The names of the list must be the names of the functional groups.  Default value is from 1 to 10)}
      #'  \item{"initBM"}{Boolean. True if using simple and bipartite SBM as initialisations. Default value  = TRUE}
      #'  \item{"maxiterVEM"}{Number of max. number of iterations in  the VEM. Default value  = 100}
      #'  \item{"maxiterVE"}{Number of max. number of iterations in  the VE. Default value  = 100}
      #'}
      optimize = function(estimOptions) {


        if (self$dependentNetwork == FALSE)
        {
          currentOptions <- list(
            verbosity     = 1,
            nbBlocksRange = lapply(1:self$nbLabels,function(l){c(1,10)}),
            nbCores       = 2,
            maxiterVE     = 100,
            maxiterVEM    = 100,
            initBM = TRUE
          )

          names(currentOptions$nbBlocksRange) <- private$dimlab
          ## Current options are default expect for those passed by the user
          currentOptions[names(estimOptions)] <- estimOptions

          # ----- formatting data for using GREMLINS
          listNetG <- lapply(private$Y, function(net) {
            if (inherits(net, "SimpleSBM_fit"))
              return(GREMLINS::defineNetwork(net$networkData,
                                    ifelse(net$directed, "diradj", "adj"),
                                    rowFG = net$dimLabels,
                                    colFG = net$dimLabels))
            if (inherits(net, "BipartiteSBM_fit"))
              return(GREMLINS::defineNetwork(net$networkData,
                                    "inc",
                                    rowFG = net$dimLabels[[1]],
                                    colFG = net$dimLabels[[2]]))
          })

          vdistrib <- private$model

          v_Kmin  <- sapply(1:self$nbLabels, function(k){currentOptions$nbBlocksRange[[k]][1]})
          v_Kmax  <- sapply(1:self$nbLabels, function(k){currentOptions$nbBlocksRange[[k]][2]})

          verbose <- (currentOptions$verbosity > 0)
          nbCores <- currentOptions$nbCores
          maxiterVE <- currentOptions$maxiterVE
          maxiterVEM <- currentOptions$maxiterVEM
          namesFG <- names(currentOptions$nbBlocksRange)
          initBM <- currentOptions$initBM

          if ( sum(abs(v_Kmin - v_Kmax)) > 0) {
            private$GREMLINSobject <- GREMLINS::multipartiteBM(
              list_Net = listNetG,
              v_distrib = vdistrib ,
              namesFG = namesFG,
              v_Kmin = v_Kmin  ,
              v_Kmax = v_Kmax ,
              v_Kinit = NULL ,
              initBM = initBM,
              keep = TRUE ,
              verbose = verbose,
              nbCores = nbCores,
              maxiterVE =  maxiterVE ,
              maxiterVEM =  maxiterVEM)
            private$import_from_GREMLINS()
          } else {
            private$GREMLINSobject <- GREMLINS::multipartiteBMFixedModel(
              list_Net = listNetG,
              v_distrib = vdistrib,
              namesFG = namesFG ,
              v_K = v_Kmax,
              classifInit = NULL,
              nbCores = nbCores,
              maxiterVE = maxiterVE,
              maxiterVEM = maxiterVEM,
              verbose = verbose)
            private$import_from_GREMLINS()
          }
        }

        else {
          currentOptions <- list(
            verbosity     = 3,
            plot          = TRUE,
            explorFactor  = 1.5,
            nbBlocksRange = c(4,Inf),
            nbCores       = 2,
            fast          = TRUE
          )
          currentOptions[names(estimOptions)] <- estimOptions

          ## Transform estimOptions to a suited for blockmodels list of options
          blockmodelsOptions <- list(
            verbosity          = currentOptions$verbosity,
            plotting           = if(currentOptions$plot) character(0) else "",
            explore_min        = currentOptions$nbBlocksRange[1],
            explore_max        = currentOptions$nbBlocksRange[2],
            ncores             = currentOptions$nbCores,
            exploration_factor = currentOptions$explorFactor
          )

          if (self$modelName[1]=="bernoulli") {model_type="bernoulli_multiplex"}
          if (self$modelName[1]=="gaussian") {model_type="gaussian_multivariate"}
          ## generating arguments for blockmodels call

          # membership type
          if (length(self$dimLabels)>1) {
            membership <-  "LBM" ; type="bipartite"
          } else {
            membership <- ifelse(!self$directed[1], "SBM_sym", "SBM") ; type="simple"
          }

          # recuperer les matrices
          Ys <- .na2zero(lapply(self$networkData,function(net) net$networkData))

          args <- list(membership_type =  membership, adj = Ys)
          args <- c(args, blockmodelsOptions)
          private$BMobject <- do.call(paste0("BM_", model_type), args)
          ## performing estimation
          private$BMobject$estimate()
          ## Exporting blockmodels output to simpleSBM_fit fields
          if (type=="simple") private$import_from_BM_Simple() else private$import_from_BM_Bipartite()

          invisible(private$BMobject)

        }
      }
  ),
  active = list(
    #' @field dependentNetwork : connection parameters in each network
    dependentNetwork = function(value) {private$dependent}
  )
)
