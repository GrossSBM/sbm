#' R6 Class definition of an Bipartite SBM fit
#'
#' This class is designed to give a representation and adjust an LBM fitted with blockmodels.
#'
#' @import R6 blockmodels
#' @export
BipartiteSBM_fit <-
  R6::R6Class(
    classname = "BipartiteSBM_fit",
    inherit = BipartiteSBM,
    private = list(
      J = NULL, # approximation of the log-likelihood
      vICL = NULL, # approximation of the ICL
      BMobject = NULL, # blockmodels output (used to stored the optimization results when blockmodels is used)
      import_from_BM = function(index = which.max(private$BMobject$ICL)) {
        private$J <- private$BMobject$PL[index]
        private$vICL <- private$BMobject$ICL[index]
        parameters <- private$BMobject$model_parameters[[index]]
        private$beta <- parameters$beta ## NULL if no covariates
        private$theta <- switch(private$BMobject$model_name,
          "bernoulli"                 = list(mean = parameters$pi),
          "bernoulli_covariates"      = list(mean = .logistic(parameters$m)),
          "bernoulli_covariates_fast" = list(mean = .logistic(parameters$m)),
          "poisson"                   = list(mean = parameters$lambda),
          "poisson_covariates"        = list(mean = parameters$lambda),
          "gaussian"                  = list(mean = parameters$mu, var = parameters$sigma2),
          "gaussian_covariates"       = list(mean = parameters$mu, var = parameters$sigma2),
          "ZIgaussian"                = list(mean = parameters$mu, var = parameters$sigma2, p0 = parameters$p0),
        )

        private$B <- vector('list',2)
        private$pi <- list(row = private$BMobject$memberships[[index]]$alpha1, col = private$BMobject$memberships[[index]]$alpha2)
        if (length(private$BMobject$memberships[[index]][["B"]]) > 0) {
          private[["B"]][[1]] <- private$BMobject$memberships[[index]][["B"]]
          if(ncol(private$Xnodes[[1]])==0){
            private[["B"]][[1]]=matrix(0,0,0)
            private$pi[[1]] = private$pi[[1]][1,]
            }
        }
        if (length(private$BMobject$memberships[[index]][["G"]]) > 0) {
          private[["B"]][[2]] <- private$BMobject$memberships[[index]][["G"]]
          if(ncol(private$Xnodes[[2]])==0){
            private[["B"]][[2]]=matrix(0,0,0)
            private$pi[[2]] = private$pi[[2]][1,]}
        }
        private$Z <- list(
          row = private$BMobject$memberships[[index]]$Z1,
          col = private$BMobject$memberships[[index]]$Z2
        )


      }
    ),
    public = list(
      #' @description constructor for a Bipartite SBM fit
      #' @param incidenceMatrix rectangular (weighted) matrix
      #' @param model character (\code{'bernoulli'}, \code{'poisson'}, \code{'gaussian'})
      #' @param dimLabels labels of each dimension (in row, in columns)
      #' @param covarList an  optional list of covariates, each of whom must have the same dimension as \code{incidenceMatrix}
      #' @param nodesCovarList an  optional list of two matrices of covariates, each of whom must have the same dimension as \code{incidenceMatrix}
      initialize = function(incidenceMatrix, model, dimLabels = c(row = "row", col = "col"), covarList = list(), nodesCovarList = vector('list',2)) {
        ## SANITY CHECKS on data
        stopifnot(is.matrix(incidenceMatrix)) # must be a matrix
        stopifnot(all(sapply(covarList, nrow) == nrow(incidenceMatrix))) # consistency of the covariates
        stopifnot(all(sapply(covarList, ncol) == ncol(incidenceMatrix))) # with the network data
        stopifnot("Nodes covariates are either not provided or a list of one or two matrices named 'row' or 'col', with as many rows as there is row or col nodes." = (length(nodesCovarList) == 0 ||
        (length(nodesCovarList) == 2  &&
        ((length(nodesCovarList[[1]]) == 0 || is.matrix(nodesCovarList[[1]]) && nrow(nodesCovarList[[1]]) == nrow(incidenceMatrix)) &&
        (length(nodesCovarList[[2]]) == 0 || is.matrix(nodesCovarList[[2]]) && nrow(nodesCovarList[[2]]) == ncol(incidenceMatrix))))))



        isBinary <- all(.na2zero(incidenceMatrix) %in% c(0, 1))
        anyRealnumber <- any(.na2zero(incidenceMatrix)%%1!=0)

        if((isBinary)&(model!="bernoulli")){stop('Choose the bernoulli distribution for your binary data')}
        if((!isBinary)&(model=="bernoulli")){stop('The bernoulli distribution is not adatped to your non binary data')}
        if(anyRealnumber&(model %in% c("poisson","bernoulli"))){stop('Choose a distribution adapted to real numbers')}



        ## INITIALIZE THE SBM OBJECT ACCORDING TO THE DATA
        connectParam <- switch(model,
          "bernoulli"  = list(mean = matrix(0, 0, 0)),
          "poisson"    = list(mean = matrix(0, 0, 0)),
          "gaussian"   = list(mean = matrix(0, 0, 0), var = 1),
          "ZIgaussian" = list(mean = matrix(0, 0, 0), var = 1, p0 = 0),
        )

        ## INITIALIZE THE SBM OBJECT ACCORDING TO THE DATA
        super$initialize(
          model = model,
          nbNodes = dim(incidenceMatrix),
          blockProp = rep(list(vector("numeric", 0)), 2),
          connectParam = connectParam,
          dimLabels = dimLabels,
          covarList = covarList,
          nodesCovarList = nodesCovarList
        )
        private$Y <- incidenceMatrix
      },
      #' @description function to perform optimization
      #' @param estimOptions a list of parameters controlling the inference algorithm and model selection. See details.
      #' @inherit estimateSimpleSBM details
      optimize = function(estimOptions = list()) {

        if (private$model == "ZIgaussian") stop("Inference not  yet  implemented for Bipartite ZI gaussian network")

        currentOptions <- list(
          verbosity = 3,
          plot = TRUE,
          exploreFactor = 1.5,
          exploreMin = 4,
          exploreMax = Inf,
          nbBlocksRange = c(4, Inf),
          nbCores = 2,
          fast = TRUE
        )
        currentOptions[names(estimOptions)] <- estimOptions

        ## Transform estimOptions to a suited for blockmodels list of options
        blockmodelsOptions <- list(
          verbosity          = currentOptions$verbosity,
          plotting           = if (currentOptions$plot) character(0) else "",
          explore_min        = currentOptions$exploreMin,
          explore_max        = currentOptions$exploreMax,
          ncores             = currentOptions$nbCores,
          exploration_factor = currentOptions$exploreFactor
        )
        fast <- currentOptions$fast

        ## generating arguments for blockmodels call

        args <- list(membership_type = "LBM", adj = private$Y)
        if (self$nbCovariates > 0) args$covariates <- private$X

        Xnodes <- private$Xnodes
        if (any(self$nbNodesCovariates > 0)){
          if (self$nbNodesCovariates[1]==0){Xnodes[[1]] <- matrix(1,private$dim[1],1)}
          if (self$nbNodesCovariates[2]==0){Xnodes[[2]] <- matrix(1,private$dim[2],1)}
          args$nodes_covariates <- Xnodes
        }


        args <- c(args, blockmodelsOptions)

        ## model construction

        model_type <- ifelse(self$nbCovariates > 0, paste0(private$model, "_covariates"), private$model)
        if (model_type == "bernoulli_covariates" & fast == TRUE) model_type <- "bernoulli_covariates_fast"
        private$BMobject <- do.call(paste0("BM_", model_type), args)

        ## performing estimation
        private$BMobject$estimate()

        ## Exporting blockmodels output to BipartiteSBM_fit fields
        private$import_from_BM()

        invisible(private$BMobject)
      },
      #' @description method to select a specific model among the ones fitted during the optimization.
      #'  Fields of the current SBM_fit will be updated accordingly.
      #' @param index integer, the index of the model to be selected (row number in storedModels)
      setModel = function(index) {
        stopifnot(!is.null(private$BMobject))
        stopifnot(index %in% self$storedModels$indexModel)
        private$import_from_BM(index)
        self$reorder()
      },
      #' @description permute group labels by order of decreasing probability
      reorder = function() {

        are_covariates <-sapply(1:2,function(l){ (self$nbNodesCovariates[l] > 0 && self$nbBlocks[l] >= 2L && is.matrix(private$pi[[l]]))})
        o  <- lapply(2:1,function(l){
         if(are_covariates[l]) {
            pi_l <- colMeans(private$pi[[l]])
          }else{
            pi_l <- private$pi[[l]]
          }
          if (l==1){
            o_l <- order( pi_l%*%private$theta$mean, decreasing = TRUE)
            }else{
            o_l <- order( private$theta$mean %*% pi_l, decreasing = TRUE)
            }
          return(o_l)})


        private$theta$mean <- private$theta$mean[o[[1]], o[[2]], drop = FALSE]
        for (l in 1:2){
          if(are_covariates[l]){
            private$pi[[l]] = private$pi[[l]][o[[l]], ]
          }else{
            private$pi[[l]] = private$pi[[l]][o[[l]]]
          }
          private$Z[[l]] <- private$Z[[l]][, o[[l]], drop = FALSE]
          if (self$nbNodesCovariates[l] > 0 && self$nbBlocks[l] >= 2L) {
            private$B[[l]] <- private$B[[l]][, o[[l]], drop = FALSE]
            private$B[[l]] <- private$B[[l]] - private$B[[l]][,ncol(private$B[[l]])]
          }
        }



        #ncol_Xnodes <- sapply(private$Xnodes,function(Mat){ifelse(is.null(Mat),0,ncol(Mat))})
        #blockProp_estim <- self$blockProp
        #nbBlocks <- self$nbBlocks
        # if(nbBlocks[1]>1){
        #   oRow <- order(private$theta$mean %*%blockProp_estim[[2]], decreasing = TRUE)
        #   }else{oRow=c(1)
        # }
        # if(nbBlocks[2]>1){
        #   oCol <- order(blockProp_estim[[1]] %*% private$theta$mean, decreasing = TRUE)
        # }else{oCol=c(1)
        # }
        #
        # if (ncol_Xnodes[1] > 0 && self$nbBlocks[1] >= 2L) {
        #     private$pi[[1]] <- private$pi[[1]][, oRow, drop = FALSE]
        #   } else {
        #     private$pi[[1]] <- private$pi[[1]][oRow]
        #   }
        #
        #
        #   if (ncol_Xnodes[2] > 0 && self$nbBlocks[2] >= 2L) {
        #     private$pi[[2]] <- private$pi[[2]][,oCol, drop = FALSE]
        #   } else {
        #   private$pi[[2]] <- private$pi[[2]][oCol]
        #   }
        #
        #
        # private$theta$mean <- private$theta$mean[oRow, oCol, drop = FALSE]
        # private$Z[[1]] <- private$Z[[1]][, oRow, drop = FALSE]
        # private$Z[[2]] <- private$Z[[2]][, oCol, drop = FALSE]
        # if (length(private$B[[1]]) > 0){
        #   private$B[[1]] <- private$B[[1]][, oRow, drop = FALSE]
        #   private$B[[1]] <- private$B[[1]] - private$B[[1]][,ncol(private$B[[1]])]
        # }
        # if (length(private$B[[2]]) > 0){
        #   private$B[[2]] <- private$B[[2]][, oCol, drop = FALSE]
        #   private$B[[2]] <- private$B[[2]] - private$B[[2]][,ncol(private$B[[2]])]
        # }
      },
      #' @description show method
      #' @param type character used to specify the type of SBM
      show = function(type = "Fit of a Bipartite Stochastic Block Model") {
        super$show(type)
        cat("* Additional fields\n")
        cat("  $probMemberships, $loglik, $ICL, $storedModels, \n")
        cat("* Additional methods \n")
        cat("  predict, fitted, $setModel, $reorder \n")
      }
    ),
    active = list(
      #' @field loglik double: approximation of the log-likelihood (variational lower bound) reached
      loglik = function(value) {
        private$J
      },
      #' @field ICL double: value of the integrated classification log-likelihood
      ICL = function(value) {
        private$vICL
      },
      #' @field penalty double, value of the penalty term in ICL
      penalty = function(value) {
        unname((self$nbConnectParam + self$nbCovariates) * log(self$nbDyads) + max(1,self$nbNodesCovariates[1])*(self$nbBlocks[1] - 1) * log(private$dim[1]) + max(1,self$nbNodesCovariates[2])*(self$nbBlocks[2] - 1) * log(private$dim[2]))
      },
      #' @field entropy double, value of the entropy due to the clustering distribution
      entropy = function(value) {
        -sum(.xlogx(private$Z[[1]])) - sum(.xlogx(private$Z[[2]]))
      },
      #' @field storedModels data.frame of all models fitted (and stored) during the optimization
      storedModels = function(value) {

        #browser()
        rowBlocks <- c(0, unlist(sapply(private$BMobject$memberships, function(m) ncol(m$Z1))))
        colBlocks <- c(0, unlist(sapply(private$BMobject$memberships, function(m) ncol(m$Z2))))
        nbConnectParam <- c(NA, unlist(sapply(private$BMobject$model_parameters, function(param) param$n_parameters)))
        U <- data.frame(
          indexModel = 1:length(private$BMobject$ICL),
          nbParams = nbConnectParam + max(1,self$nbNodesCovariates[1])*(rowBlocks-1) + max(1,self$nbNodesCovariates[2])*(colBlocks - 1),
          nbRowBlocks = rowBlocks,
          nbColBlocks = colBlocks,
          ICL = c(private$BMobject$ICL),
          loglik = c(private$BMobject$PL)
        )
        U[!is.na(U$nbParams), , drop = FALSE]
      }
    )
  )
