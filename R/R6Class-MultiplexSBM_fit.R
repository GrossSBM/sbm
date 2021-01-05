#' R6 Class definition of a Multiplex SBM fit
#'
#' This class is designed to give a representation and adjust a Multiplex SBM fitted with GREMLIN.
#'
#' @import R6 GREMLINS
#' @export
MultiplexSBM_fit <-
  R6::R6Class(
    classname = "MultiplexSBM_fit",
    inherit = MultiplexSBM,
    # fields for internal use (referring to the mathematical notation)
    private = list(
      BMobject = NULL,
      GREMLINSobject       = NULL,

      #------------ function to convert GREMLINS result into a sbm object result
      import_from_GREMLINS = function(index = 1) {

        GREMLINSfit <- private$GREMLINSobject$fittedModel[[index]]
        list_pi <- lapply(private$namesFG,function(n_){GREMLINSfit$paramEstim$list_pi[[n_]]})
        list_tau <- lapply(private$namesFG,function(n_){GREMLINSfit$paramEstim$tau[[n_]]})
        list_theta <-lapply(1:private$nbNet, function(s_){
          GREMLINSfit$paramEstim$list_theta[[paste(private$namesFG[private$E[s_,1]],private$namesFG[private$E[s_,2]],sep='')]]
        })
        #-----------------------------------------------------
        list_theta_mean <- lapply(1:private$nbNet,function(s_){
          if(private$distrib[s_] %in% c('ZIgaussian','gaussian')){u = list_theta[[s_]]$mean}else{u = list_theta[[s_]]}
          u})

        ordAll <- order_mbm(list_theta_mean,list_pi,private$E)
        #----------------------------------------------------------
        listtau <- lapply(1:private$nbFG, FUN = function(s){list_tau[[s]][,ordAll[[s]]]})
        listpi <-  lapply(1:private$nbFG, FUN = function(s){list_pi[[s]][ordAll[[s]]]})
        names(listpi)<- names(listtau)<- private$namesFG

        lapply(private$listNet, function(net) {
          if (substr(class(net)[1], 1, 6) == "Simple") {
            Lab <- net$dimLabels[[1]]
            net$varProb <- listtau[[Lab]]
            net$blockProp <- listpi[[Lab]]
          }
          else {
            rowLab <- net$dimLabels[[1]]
            colLab <- net$dimLabels[[2]]
            net$varProb <-
              list(listtau[[rowLab]], listtau[[colLab]])
            net$blockProp <-
              list(listpi[[rowLab]], listpi[[colLab]])
          }
        })

        #print(private$nbNet)
        private$theta = list()
        for (s_ in 1:private$nbNet){
          o_row <- ordAll[[private$E[s_,1]]]
          o_col <- ordAll[[private$E[s_,2]]]
          l_s <- list(mean  =  list_theta_mean[[s_]][o_row ,o_col])
          if (private$distrib[s_] %in% c('gaussian','ZIgaussian')){
            var_s <- list_theta[[s_]]$var
            if (is.matrix(var_s)){var_s  = var_s[o_row,o_col]}
            l_s$var <- var_s
          }
          if (private$distrib[s_] == 'ZIgaussian'){
            p0_s <- list_theta[[s_]]$p0
            if (is.matrix(p0_s)){p0_s  = p0_s[o_row,o_col]}
            l_s$p0 <- p0_s
          }


          private$theta[[s_]]=l_s
          private$listNet[[s_]]$connectParam = l_s
        }

        private$tau = listtau

        private$allZ = lapply(1:length(listtau),function(l){
          if(!is.matrix(listtau[[l]])){listtau[[l]] = matrix(listtau[[l]],ncol=1)}
          u <- apply(listtau[[l]],1,which.max)
          u})
        private$pi  = listpi
      },
      import_from_BM_Simple = function(index = which.max(private$BMobject$ICL)) { # a function updating the Class
        super$import_from_BM(index)
        private$tau <- private$BMobject$memberships[[index]]$Z
        private$pi  <- colMeans(private$tau)
      },
      import_from_BM_Bipartite  = function(index = which.max(private$BMobject$ICL)) {
        super$import_from_BM(index)
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
      #' @param listSBM list of SBM object with
      #' @param dep boolean indicating whether dependence is assumed between networks beyond the common dependence on the latent variables
      initialize = function(listSBM,dep=FALSE) {
        super$initialize(listSBM,dep)
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


        if (self$modelDependence==FALSE)
        {
        currentOptions <- list(
          verbosity     = 1,
          nbBlocksRange = lapply(1:private$nbFG,function(l){c(1,10)}),
          nbCores       = 2,
          maxiterVE     = 100,
          maxiterVEM    = 100,
          initBM = TRUE
        )
        names(currentOptions$nbBlocksRange) <- private$namesFG
        ## Current options are default expect for those passed by the user
        currentOptions[names(estimOptions)] <- estimOptions



        # ----- formatting data for using GREMLINS
        listNetG <- lapply(private$listNet, function(net) {
          if (substr(class(net)[1], 1, 6) == "Simple") {
            ifelse(net$directed, type <- "diradj", type <- "adj")
          }
          else {
            type <-  "inc"
          }
          GREMLINS::defineNetwork(net$netMatrix,
                                  type,
                                  rowFG = net$dimLabels[[1]],
                                  colFG = net$dimLabels[[2]])
        })


        vdistrib <- private$distrib
        v_Kmin  <- sapply(1:private$nbFG, function(k){currentOptions$nbBlocksRange[[k]][1]})
        v_Kmax  <- sapply(1:private$nbFG, function(k){currentOptions$nbBlocksRange[[k]][2]})
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
          fast  = currentOptions$fast

          if (self$modelName[1]=="bernoulli") {model_type="bernoulli_multiplex"}
          if (self$modelName[1]=="bernoulli") {model_type="gaussian_multivariate"}
          ## generating arguments for blockmodels call

          #membershiptype =

          # recuperer les matrices

          args <- list(membership_type =  ifelse(!private$directed_, "SBM_sym", "SBM"), adj = .na2zero(private$Y))
          args <- c(args, blockmodelsOptions)


          private$BMobject <- do.call(paste0("BM_", model_type), args)

        }


      }

    )
  )
