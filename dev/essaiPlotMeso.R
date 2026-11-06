library(sbm)
data("multipartiteEcologicalNetwork")
load('/home/sophie/WORK_LOCAL/RECHERCHE/FORMATIONS_RECHERCHE/FORMATIONS_RESEAUX/FORMATION_CESAB_ECONET/econetoolbox/data/day4_files/res_Multipartite_BM_Dattilo.rda')

library(RColorBrewer)
res <- res_MSBM_dattilo

nbNet = length(res$connectParam)
list_pi <-res$blockProp
theta = lapply(1:nbNet,function(l){res$connectParam[[l]]$mean})

E <- res$architecture

FG <- substr(res$dimLabels,1,3)
nFG = length(FG)
nbBlocks <- res$nbBlocks




meso_netDF <-do.call(rbind, lapply(1:nbNet,function(n){
  E1 <- E[n,1]
  E2 <- E[n,2]
  FG1 <- FG[E1]
  FG2 <- FG[E2]
  nB1 <- nbBlocks[E1]
  nB2 <- nbBlocks[E2]
  theta_n <- theta[[n]]
  theta_n <- theta_n#*(theta_n>0.05)
  row.names(theta_n) <- paste(rep(FG1,nB1),1:nB1,sep=' ')
  colnames(theta_n) <- paste(rep(FG2,nB2),1:nB2,sep=' ')
  DF_n <- as.data.frame(as.table(theta_n),stringsAsFactors = FALSE)
  names(DF_n) <- c('from','to','weight')
  return(DF_n)
}))



meso_net <- as.network(meso_netDF)
meso_net  %v% "FG" <- as.character(sub(" .*", "", network.vertex.names(meso_net)))



mypalette <- brewer.pal(n = length(unique(meso_net %v% "FG")), name = "Set2")
fg_levels <- unique(meso_net %v% "FG")
colors <- setNames(brewer.pal(length(fg_levels), "Set1"), fg_levels)

meso_net %v% "blocksize" <- unlist(lapply(1:nFG,function(k){res$blockProp[[k]]*res$nbNodes[k]}))
meso_net %e% "plotsize" <- 2*meso_net %e% "weight"
meso_net %v% "plotNodesize" <- 1000*meso_net %v% "blocksize"

ggnet2(
  meso_net,
  color = "FG",
  palette = colors,
  size = "plotNodesize",
  label = TRUE,         # show node labels
  label.color = "white",# label color
  label.size = 3,
  edge.size = "plotsize",
) + guides(size = "none", color= "none")


