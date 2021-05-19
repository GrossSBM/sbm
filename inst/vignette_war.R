rm(list=ls())
library(missSBM)
data(war)
library(igraph)

V(war$alliance)%in%V(war$belligerent)


A = as.matrix(get.adjacency(war$alliance))
B = matrix(0,nrow(A),ncol(A))
B[1:length(V(war$belligerent)),1:length(V(war$belligerent))] = as.matrix(get.adjacency(war$belligerent))
plot(graph_from_adjacency_matrix(B))
sum(rowSums(B)>0)
sum(diag(B))
isSymmetric(B)
# sous partie
B=B[1:83,1:83]
A=A[1:83,1:83]

library(sbm)
netA = defineSBM(A,model="bernoulli",dimLabels = "country") # bug dimLabels par defaut
netB = defineSBM(B,model="bernoulli",dimLabels = "country")
plotMyMultiplexMatrix(list(netA,netB))


MultiplexFitdep = estimateMultiplexSBM(list(netA,netB),dependent = T)
plot(MultiplexFitdep,type="data")
plot(MultiplexFitdep,type="expected")

MultiplexFitdep$memberships
p11 = MultiplexFitdep$connectParam$prob11
p01 = MultiplexFitdep$connectParam$prob01
p10 = MultiplexFitdep$connectParam$prob10
# proba de faire la guerre sachant alliance
round(p11/(p11+p10),1)
# proba de faire la guerre
round(p11+p01,1)

round(p11/(p01+p11),1)
image(p01+p11)
image(p10+p11)
image(p11/(p11+p10))

MultiplexFitIndep = estimateMultiplexSBM(list(netA,netB),dependent = F)
plot(MultiplexFitIndep,type="expected")

MultiplexFitdep$ICL
MultiplexFitIndep$ICL
#save.image("inst/resWARallnetwork.Rdata")

