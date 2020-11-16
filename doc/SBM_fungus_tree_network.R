## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, message=FALSE, warning=FALSE--------------------------------------
library(sbm)
library(ggplot2)
library(knitr)
theme_set(theme_bw())

## ----import dataset-----------------------------------------------------------
data("fungusTreeNetwork")
str(fungusTreeNetwork, max.level = 1)

## ----tree_tree_binary network-------------------------------------------------
tree_tree_binary <- 1 * (fungusTreeNetwork$tree_tree != 0)

## ----tree_tree_binary network plot data---------------------------------------
plotMyMatrix(tree_tree_binary, dimLabels = list(row = 'tree', col = 'tree') )

## ----simpleSBM----------------------------------------------------------------
mySimpleSBM <- tree_tree_binary %>% 
  estimateSimpleSBM("bernoulli", dimLabels = list(row = 'tree', col='tree'), estimOptions = list(verbosity = 0, plot = FALSE))

## ----simpleSBMfit-------------------------------------------------------------
class(mySimpleSBM)
mySimpleSBM

## ----impleSBMfit fields-------------------------------------------------------
mySimpleSBM$nbBlocks
mySimpleSBM$nbNodes
mySimpleSBM$nbCovariates

## ----simpleSBMfit plot1-------------------------------------------------------
plot(mySimpleSBM, type = "data", dimLabels  = list(row = 'tree', col= 'tree'))

## ----simpleSBMfit plot2-------------------------------------------------------
plot(mySimpleSBM, type = "expected")

## ----simpleSBMfit plotmeso----------------------------------------------------
plot(mySimpleSBM, type = "meso")

## ----simpleSBMfit coef--------------------------------------------------------
coef(mySimpleSBM, 'block')
coef(mySimpleSBM, 'connectivity')

## ----simpleSBM storedModel----------------------------------------------------
mySimpleSBM$storedModels %>% kable()

## ----simpleSBM ICL------------------------------------------------------------
mySimpleSBM$storedModels %>%  ggplot() + aes(x = nbBlocks, y = ICL)  + geom_line() + geom_point(alpha = 0.5)

## ----simpleSBMfit changeModel-------------------------------------------------
mySimpleSBM$setModel(4)
mySimpleSBM$nbBlocks
mySimpleSBM$plot(type = 'expected')
mySimpleSBM$setModel(5)

## ----tree_tree network plot data----------------------------------------------
tree_tree <- fungusTreeNetwork$tree_tree
plotMyMatrix(tree_tree, dimLabels = list(row = 'tree', col = 'tree'))

## ----simpleSBM Poisson--------------------------------------------------------
mySimpleSBMPoisson <- tree_tree %>% 
  estimateSimpleSBM("poisson", dimLabels = list(row = 'tree', col= 'tree'),estimOptions = list(verbosity = 0, plot = FALSE))

## ----simpleSBMfitPoisson------------------------------------------------------
class(mySimpleSBMPoisson)
mySimpleSBMPoisson

## ----impleSBMfitPoison fields-------------------------------------------------
mySimpleSBMPoisson$nbBlocks
mySimpleSBMPoisson$nbNodes
mySimpleSBMPoisson$nbCovariates

## ----simpleSBMfitPoisson plot1------------------------------------------------
plot(mySimpleSBMPoisson, type = "data", dimLabels = list(row = 'tree', col= 'tree'))

## ----simpleSBMfitPoisson plot2------------------------------------------------
plot(mySimpleSBMPoisson, type = "expected", dimLabels = list(row = 'tree', col= 'tree'))

## ----simpleSBMfitPoisson plotmeso---------------------------------------------
plot(mySimpleSBMPoisson, type = "meso")

## ----simpleSBMfitPoisson coef-------------------------------------------------
coef(mySimpleSBMPoisson, 'block')
coef(mySimpleSBMPoisson, 'connectivity')

## ----covar SBM,echo=TRUE,eval= TRUE-------------------------------------------
mySimpleSBMCov<- 
  tree_tree %>% 
  estimateSimpleSBM(
    model = 'poisson', 
    directed = FALSE, 
    dimLabels = list(row = 'tree', col= 'tree'),
    covariates  = fungusTreeNetwork$covar_tree, 
    estimOptions = list(verbosity = 0, plot = FALSE, nbCores = 2)
  )

## ----select SBM covar, echo=TRUE, eval = TRUE---------------------------------
mySimpleSBMCov$nbBlocks

## ----extract param SBM poisson covar, echo=TRUE, eval = TRUE------------------
mySimpleSBMCov$connnectParam
mySimpleSBMCov$blockProp
mySimpleSBMCov$memberships
mySimpleSBMCov$covarParam

## ----simpleSBMfitPoisson covar coef-------------------------------------------
coef(mySimpleSBMCov, 'covariates')

## ----simpleSBMfitPoisson covar fitted, results='hide'-------------------------
fitted(mySimpleSBMCov)
predict(mySimpleSBMCov)
predict(mySimpleSBMCov, fungusTreeNetwork$covar_tree)

## ----plot incidence-----------------------------------------------------------
plotMyMatrix(fungusTreeNetwork$fungus_tree, dimLabels = list(row = 'fungis', col= 'tree'))

## ----tree_fungi_bipartite network---------------------------------------------
myBipartiteSBM <- 
  fungusTreeNetwork$fungus_tree %>% 
  estimateBipartiteSBM(model = 'bernoulli', dimLabels = list(row = 'fungis', col= 'tree'),estimOptions = list(verbosity = 0, plot = FALSE))

## ----bipartite.sbm fields-----------------------------------------------------
myBipartiteSBM$nbNodes
myBipartiteSBM$nbBlocks
myBipartiteSBM$connectParam
coef(myBipartiteSBM, 'block')
coef(myBipartiteSBM, 'connectivity')

## ----plot bipartite-----------------------------------------------------------
plot(myBipartiteSBM, dimLabels = list(row = 'fungis', col = 'tree'))

## ----plot bipartite meso------------------------------------------------------
plot(myBipartiteSBM, type  = 'meso', 
     dimLabels = list(row = 'fungis', col= 'tree'),
     plotOptions = list(edge.width = 1, vertex.size = c(1,2)))

