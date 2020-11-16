## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
collapse = TRUE,
comment = "#>"
)

## ----setup, message=FALSE, warning=FALSE--------------------------------------
library(sbm)
library(ggplot2)
library(knitr)

## ----loading dataset, eval=TRUE-----------------------------------------------
data(multipartiteEcologicalNetwork)
str(multipartiteEcologicalNetwork)
names(multipartiteEcologicalNetwork)

## ----transform dataset,  eval=TRUE--------------------------------------------
Net <- multipartiteEcologicalNetwork
type='bipartite'
model = 'bernoulli'
directed = FALSE
PlantFlovis = defineSBM(Net$Inc_plant_flovis, model,type,directed,
                        dimLabels = list(row="Plants",col="Flovis"))
PlantAnt = defineSBM(Net$Inc_plant_ant,model,type,directed,
                      dimLabels =list(row = "Plants", col = "Ants"))
PlantBird = defineSBM(Net$Inc_plant_bird,model,type,directed,
                     dimLabels =list(row = "Plants",col = "Birds"))

## ----example of dataset, eval=TRUE--------------------------------------------
PlantFlovis$netMatrix[1:2,1:2]

## ----plot data----------------------------------------------------------------
plotMyMultipartiteMatrix(list(PlantFlovis,PlantAnt,PlantBird))

## ----load result, echo = FALSE, eval = TRUE-----------------------------------
load('resMultipartiteEcological.rda')

## ----MBM, echo = TRUE, eval = FALSE-------------------------------------------
#  estimOptions = list(initBM = FALSE)
#  listSBM <- list(PlantFlovis, PlantAnt, PlantBird)
#  myMSBM <- estimateMultipartiteSBM(listSBM,estimOptions)

## ----MBM what-----------------------------------------------------------------
myMSBM

## ----MBM v_K------------------------------------------------------------------
myMSBM$nbBlocks

## ----MBM param----------------------------------------------------------------
myMSBM$blockProp
myMSBM$connectParam

## ----MBM Z--------------------------------------------------------------------
table(myMSBM$memberships$Plants)
table(myMSBM$memberships$Ants)      

## ----storedmodels-------------------------------------------------------------
myMSBM$storedModels

## ----plot, eval = TRUE--------------------------------------------------------
plot(myMSBM) 

## ----plot meso, eval = FALSE--------------------------------------------------
#  plot(myMSBM,type = 'meso',plotOptions=list(vertex.size = 0.5))

