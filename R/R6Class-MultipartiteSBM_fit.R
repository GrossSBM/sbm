#' R6 Class definition of a Multipartite SBM fit
#'
#' This class is designed to give a representation and adjust a Multipartite SBM fitted with GREMLIN.
#'
#' @import R6 GREMLINS
#' @export
MultipartiteSBM_fit <-
  R6::R6Class(
    classname = "MultipartiteSBM_fit",
    inherit = MultipartiteSBM,
    # fields for internal use (referring to the mathematical notation)
    private = list(
       tau            = NULL, # variational parameters for posterior probability of class belonging
       GREMLINSobject = NULL,

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
      }
    ),
    #-----------------------------------------------
    public = list(
      #' @description constructor for Multipartite SBM
      #' @param netList list of SBM objects
      initialize = function(netList) {
        private$netList <- netList
        model     <- map_chr(listSBM, "modelName")
        directed  <- map(netList, "directed") %>% map_lgl(~ifelse(is.null(.x), NA, .x))
        dimension <- netList %>% map("dimension") %>% unlist() %>% unique()
        dimLabels <- netList %>% map("dimLabels") %>% unlist() %>% unique()

        arch      <- listSBM %>% map_df("dimLabels") %>%
           map(factor, levels = private$namesFG) %>% map_df(as.numeric) %>% as.matrix()

        ## alternative to code above for architecture and dimension
        # ###
        # E_FG <- lapply(listSBM,function(net){return(c(net$dimLabels$row,net$dimLabels$col))})
        # E_FG <- do.call(rbind,E_FG)
        # E <- matrix(sapply(E_FG,function(a){which(private$namesFG == a)}), self$nbNetworks,2)
        # private$E <-  E
        # private$dimFG <- sapply(1:self$nbLabels ,function(k){
        #   u <- which(E[,1] == k); v = 1;
        #   if (length(u) == 0) {u <- which(E[,2] == k); v = 2}
        #   u <- u[1]
        #   dim(listSBM[[u]]$netMatrix)[v]}
        # )
###      private$allZ <- memberships

        super$initialize(model        = model,
                         architecture = arch,
                         directed     = directed,
                         dimension    = dimension,
                         dimLabels    = dimLabels,
                         blockProp    = blockProp,
                         connectParam = connectParam)
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
            type <- ifelse(net$directed, "diradj", "adj")
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
      },
      #' @description prediction under the currently estimated model
      #' @return a list of matrices matrix of expected values for each dyad
      predict = function() {
        map(private$listNet, predict)
      },
      #' @description method to select a specific model among the ones fitted during the optimization.
      #'  Fields of the current MultipartiteSBM_fit will be updated accordingly.
      #' @param index integer, the index of the model to be selected (row number in storedModels)
      setModel = function(index) {
        stopifnot(!is.null(private$GREMLINSobject))
        stopifnot(index %in% seq.int(nrow(self$storedModels)))
        private$import_from_GREMLINS(index)
      },
      #' @description print method
      #' @param type character to tune the displayed name
      show = function(type = "Fit of a Multipartite Stochastic Block Model"){
        cat(type, "\n")
        cat(self$nbLabels, "functional groups (", self$dimLabels, "), ", self$nbNetworks, "networks\n")
        cat("=====================================================================\n")
        cat("nbNodes per FG = (", self$nbNodes, ") --  nbBlocks per FG = (",self$nbBlocks, ")\n")
        cat("distributions on each network =(", self$modelName ,")\n")
        cat("=====================================================================\n")
        cat("* Useful fields \n")
        cat(" $nbNetwork, $nbNodes, $nbBlocks, $dimLabels, $archiMultipartite \n")
        cat(" $modelName, $blockProp, $connectParam, $memberships, $probMemberships\n")
        cat("* Useful functions \n")
        cat("$plot, $optimize, $predict, $setModel, $storedModels \n")
      },
      #' @description print method
      print = function() self$show()
  ),
  #-----------------------------------------------
  active = list(
    #' @field memberships a list with the memberships in all the functional groups
    memberships = function(value) {setNames(lapply(private$tau, as_clustering), private$dimlab)}
    #' @field probMemberships or list of nbFG matrices for of estimated probabilities for block memberships for all nodes
    probMemberships = function(value) {private$tau},
    #' @field nbBlocks : vector with the number of blocks in each FG
    nbBlocks = function(value) {setNames(sapply(private$tau, ncol), private$dimlab)},
    #' @field storedModels data.frame of all models fitted (and stored) during the optimization
    storedModels = function(value) {
      GO <- private$GREMLINSobject
      nbModels <- length(GO$fittedModel)
      Blocks <- as.data.frame(t(sapply(GO$fittedModel, function(m) m$paramEstim$v_K)))
      colnames(Blocks) <- paste('nbBlocks',private$dimlab)
      nbConnectParam <- sapply(GO$fittedModel, function(m){
        computeNbConnectParams_MBM(m$paramEstim$v_K, private$model, private$arch, private$directed_)
      })
      nbParams  <- nbConnectParam + rowSums(Blocks) - ncol(Blocks)
      indexModel <- 1:nbModels
      U <- cbind(indexModel,nbParams, Blocks)
      U$nbBlocks <-rowSums(Blocks)
      U$ICL <- sapply(GO$fittedModel, function(m) m$ICL)
      U$loglik  <- sapply(GO$fittedModel,function(m){len <- length(m$vJ); m$vJ[len]})
      U
    }
  )
)

