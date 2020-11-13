#' R6 Class definition of a Multipartite SBM fit
#'
#' This class is designed to give a representation and adjust a Multipartite SBM fitted with GREMLIN.
#'
#' @import R6 GREMLIN
#' @export
MultipartiteSBM_fit <-
  R6::R6Class(
    classname = "MultipartiteSBM_fit",
    inherit = MultipartiteSBM,
    # fields for internal use (referring to the mathematical notation)
    private = list(
      GREMLINobject       = NULL,
      import_from_GREMLIN = function(index = 1) {

        GREMLINfit <- private$GREMLINobject$fittedModel[[index]]
        list_pi <- lapply(private$namesFG,function(n_){GREMLINfit$paramEstim$list_pi[[n_]]})
        list_tau <- lapply(private$namesFG,function(n_){GREMLINfit$paramEstim$tau[[n_]]})
        list_theta <-lapply(1:private$nbNet, function(s_){
          GREMLINfit$paramEstim$list_theta[[paste(private$namesFG[private$E[s_,1]],private$namesFG[private$E[s_,2]],sep='')]]
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
        }
    ),
    #-----------------------------------------------
    public = list(
      #' @description constructor for Multipartite SBM
      #' @param listSBM list of SBM object with
      initialize = function(listSBM) {
        super$initialize(listSBM)
      },
      #' @description estimation via GREMLIN
      #' @param currentOptions options for MultipartiteBM
      optimize = function(currentOptions) {


        # ----- formatting data for using GREMLIN
        listNetG <- lapply(private$listNet, function(net) {
          if (substr(class(net)[1], 1, 6) == "Simple") {
            ifelse(net$directed, type <- "diradj", type <- "adj")
          }
          else {
            type <-  "inc"
          }
          GREMLIN::defineNetwork(net$netMatrix,
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
          private$GREMLINobject <- GREMLIN::multipartiteBM(
          list_Net = listNetG,
          v_distrib = vdistrib ,
          namesFG = namesFG,
          v_Kmin = v_Kmin  ,
          v_Kmax = v_Kmax ,
          v_Kinit = NULL ,
          initBM = initBM,
          save = TRUE ,
          verbose = verbose,
          nbCores = nbCores,
          maxiterVE =  maxiterVE ,
          maxiterVEM =  maxiterVEM)
          private$import_from_GREMLIN()
        } else {
          private$GREMLINobject <- GREMLIN::multipartiteBMFixedModel(
            list_Net = listNetG,
            v_distrib = vdistrib,
            namesFG = namesFG ,
            v_K = v_Kmax,
            classifInit = NULL,
            nbCores = nbCores,
            maxiterVE = maxiterVE,
            maxiterVEM = maxiterVEM,
            verbose = verbose)
        }
        ### TODO find what to ouput???
        #import_from_GREMLIN doit remplir listFit
      },
      #' @description getBM returns a given network from a Multipartite SBM
      #' @param i index of the asked network
      getBM = function(i) {
        private$listNet[[i]]
      },
      #' @description prediction under the currently estimated model
      #' @param covarList a list of covariates. By default, we use the covariates with which the model was estimated
      #' @return a list of matrices matrix of expected values for each dyad
      predict = function() {
          lapply(1:private$nbNet,function(l){BMl <- self$getBM(l); BMl$predict()})
      },
      #' @description method to select a specific model among the ones fitted during the optimization.
      #'  Fields of the current MultipartiteSBM_fit will be updated accordingly.
      #' @param index integer, the index of the model to be selected (row number in storedModels)
      setModel = function(index) {
        stopifnot(!is.null(private$GREMLINobject))
        stopifnot(index %in% seq.int(nrow(self$storedModels)))
        private$import_from_GREMLIN(index)
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
        cat("$plot, $optimize, $predict, $setModel, $getBM, $storedModels \n")
      },
      #' @description print method
      print = function() self$show()
    ),
    #-----------------------------------------------
    active=list(
      #' @field memberships a list with the memberships in all the functional groups
    memberships = function(value) {
      if (missing(value)){
        M <- private$allZ
        names(M) <- private$namesFG
        return(M)}  else {private$allZ <- value}},
    #' @field probMemberships  or list of nbFG matrices for of estimated probabilities for block memberships for all nodes
    probMemberships = function(value) {private$tau  },
    #' @field nbBlocks : vector with the number of blocks in each FG
    nbBlocks = function(value) {
      r<- sapply(private$allZ, function(z){length(unique(z))})
      names(r) <- private$namesFG
      return(r)},
    #' @field storedModels data.frame of all models fitted (and stored) during the optimization
    storedModels = function(value) {
      GO <- private$GREMLINobject
      nbModels <- length(GO$fittedModel)
      Blocks <- as.data.frame(t(sapply(GO$fittedModel, function(m) m$paramEstim$v_K)))
      colnames(Blocks) <- paste('nbBlocks',private$namesFG)
      nbConnectParam <-sapply(GO$fittedModel, function(m){
        E <- private$E;
        distrib <- private$distrib
        directed <- private$directed_
        r <- computeNbConnectParams_MBM(m$paramEstim$v_K,distrib,E,directed)
        r})
      nbParams  <- nbConnectParam + rowSums(Blocks) - ncol(Blocks)
      indexModel <- 1:nbModels
      U <- cbind(indexModel,nbParams, Blocks)
      U$nbBlocks <-rowSums(Blocks)
      U$ICL <- sapply(GO$fittedModel, function(m) m$ICL)
      U$loglik  <- sapply(GO$fittedModel,function(m){len <- length(m$vJ); m$vJ[len]})
      return(U)
    },
    #' @field blockProp : block proportions in each function group
    blockProp = function(value) {private$pi},
    #' @field connectParam : connection parameters in each network
    connectParam = function(value) {private$theta}
)
)




