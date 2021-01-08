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
      J = NULL,
      vICL = NULL,
      BMobject = NULL,
      GREMLINSobject       = NULL,

      #------------ function to convert GREMLINS result into a sbm object result
      import_from_GREMLINS = function(index = 1) {

        GREMLINSfit <- private$GREMLINSobject$fittedModel[[index]]
        list_pi <- lapply(private$namesFG,function(n_){GREMLINSfit$paramEstim$list_pi[[n_]]})
        list_tau <- lapply(private$namesFG,function(n_){GREMLINSfit$paramEstim$tau[[n_]]})
        list_theta <-lapply(1:self$nbNetworks, function(s_){
          GREMLINSfit$paramEstim$list_theta[[paste(private$namesFG[private$E[s_,1]],private$namesFG[private$E[s_,2]],sep='')]]
        })
        #-----------------------------------------------------
        list_theta_mean <- lapply(1:self$nbNetworks,function(s_){
          if(private$model[s_] %in% c('ZIgaussian','gaussian')){u = list_theta[[s_]]$mean}else{u = list_theta[[s_]]}
          u})

        ordAll <- order_mbm(list_theta_mean,list_pi,private$E)
        #----------------------------------------------------------
        listtau <- lapply(1:self$nbLabels, FUN = function(s){list_tau[[s]][,ordAll[[s]]]})
        listpi <-  lapply(1:self$nbLabels, FUN = function(s){list_pi[[s]][ordAll[[s]]]})
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
        for (s_ in 1:self$nbNetworks){
          o_row <- ordAll[[private$E[s_,1]]]
          o_col <- ordAll[[private$E[s_,2]]]
          l_s <- list(mean  =  list_theta_mean[[s_]][o_row ,o_col])
          if (private$model[s_] %in% c('gaussian','ZIgaussian')){
            var_s <- list_theta[[s_]]$var
            if (is.matrix(var_s)){var_s  = var_s[o_row,o_col]}
            l_s$var <- var_s
          }
          if (private$model[s_] == 'ZIgaussian'){
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

      import_from_BM  = function(index = which.max(private$BMobject$ICL)) {
        private$J     <- private$BMobject$PL[index]
        private$vICL  <- private$BMobject$ICL[index]
        parameters    <- private$BMobject$model_parameters[[index]]
        private$theta <- switch(private$BMobject$model_name,
                                "gaussian_multivariate" = list(mean=parameters$mu,cov=parameters$Sigma),
                                "bernoulli_multiplex" = list(prob00=parameters$pi$`00`,prob01=parameters$pi$`01`,
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
      #' @param listSBM list of SBM object with
      #' @param dep boolean indicating whether dependence is assumed between networks beyond the common dependence on the latent variables
      initialize = function(listSBM,dep=FALSE) {
        super$initialize(listSBM,dep=dep)
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
          nbBlocksRange = lapply(1:self$nbLabels,function(l){c(1,10)}),
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
          fast  = currentOptions$fast

          if (self$modelName[1]=="bernoulli") {model_type="bernoulli_multiplex"}
          if (self$modelName[1]=="gaussian") {model_type="gaussian_multivariate"}
          ## generating arguments for blockmodels call


          # membership type
          if (length(unique(self$dimLabels))>1) {membership <-  "LBM" ; type="bipartite"}
          else {membership <- ifelse(!self$directed[1], "SBM_sym", "SBM") ; type="simple"}

          # recuperer les matrices
          Ys <- .na2zero(lapply(self$listSBM,function(net) net$netMatrix))

          args <- list(membership_type =  membership, adj = Ys)
          args <- c(args, blockmodelsOptions)
          private$BMobject <- do.call(paste0("BM_", model_type), args)
          ## performing estimation
          private$BMobject$estimate()
          print(private$BMobject)
          ## Exporting blockmodels output to simpleSBM_fit fields
          if (type=="simple") private$import_from_BM_Simple() else private$import_from_BM_Bipartite()

          invisible(private$BMobject)

        }
      }
    ),
    active = list(
      #' @field memberships a list with the memberships in all the functional groups
      memberships = function(value) {
        if (missing(value)) {
          return(setNames(private$allZ, private$namesFG))
        } else {private$allZ <- value}},
      #' @field probMemberships or list of nbFG matrices for of estimated probabilities for block memberships for all nodes
      probMemberships = function(value) {private$tau},
      #' @field nbBlocks : vector with the number of blocks in each FG
      nbBlocks = function(value) {
        setNames(sapply(private$allZ, function(z){length(unique(z))}), private$namesFG)
      },
      #' @field blockProp : block proportions in each function group
      blockProp = function(value) {private$pi},
      #' @field connectParam : connection parameters in each network
      connectParam = function(value) {private$theta}
    )
  )
