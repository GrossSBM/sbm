plotmyMultipartiteMatrix = function(list_Net){
  
  
  
  mydataR6 <- GREMLIN:::formattingData(list_Net)
  
  E <- mydataR6$E
  Q <- mydataR6$Q
  v_NQ <- mydataR6$v_NQ
  G <- matrix(0,Q,Q)
  G[mydataR6$E[,1],mydataR6$E[,2] ] = 1
  uRow <- which(rowSums(G)>0)
  uCol <- which(colSums(G)>0)
  
  cardE <- mydataR6$cardE
   
  EndBlocks <- cumsum(v_NQ)
  BegBlocks <- c(0,cumsum(v_NQ)[-Q]) +1
  
  N <- sum(v_NQ)
  MetaMat <- matrix(NA,N,N)
  dim(MetaMat)

  FGCol <- FGRow <- rep(mydataR6$namesFG, times = v_NQ)
  indFGCol <- indFGRow <- rep(1:Q, times = v_NQ)
  for (i in 1:Q){
    for (j in 1:Q){
      if (G[i,j] == 1){
      place_row <-  BegBlocks[i]: EndBlocks[i]
      place_col <-  BegBlocks[j]: EndBlocks[j]
      kij <- which(rowSums(E == matrix(rep(c(i,j),cardE),ncol=2,nrow= cardE,byrow = T)) == 2)
      
      MetaMat[place_row, place_col] <- mydataR6$mats[[kij]]
      }  
    }
  }
  uRow <- which(apply(MetaMat , 1,function(u){all(is.na(u) )}))
  uCol <- which(apply(MetaMat , 2,function(u){all(is.na(u))}))
  MetaMatSmall <- MetaMat[-uRow, -uCol]
  FGColSmall <- FGCol[-uCol]
  FGRowSmall <- FGRow[-uRow]
  indFGColSmall <- indFGCol[-uCol]
  indFGRowSmall <- indFGRow[-uRow]
  
  g <- plotMyMatrix(MetaMatSmall)  
  resPlot <- sbm:::plotMatrix(MetaMatSmall,dimLabels = NULL, clustering = list(row=indFGRowSmall, col=indFGColSmall))
  dataPlot <- resPlot$data
  dataPlot$FG_row <- as.factor(FGRowSmall [N - dataPlot$index_row])
  dataPlot$FG_col <- as.factor(FGColSmall [dataPlot$index_col])
  g <- ggplot(data = dataPlot, aes(y = index_row, x = index_col, fill = link))
  g <- g + geom_tile()
  if (!binary) {g <-  g +  scale_fill_gradient(low = "white", high = "black", limits = u,na.value = "transparent")}
  if (binary) {g <- g + scale_fill_manual(breaks = c("0", "1"),values = c("white", "black"),na.value = "transparent")}
  g <- g  +  scale_x_discrete(drop = FALSE) + scale_y_discrete(drop = FALSE)
  g <- g + theme(axis.text.x = element_text(angle = 270, hjust = 0))
  g <- g +  labs(x = colFG, y = rowFG) +  theme(aspect.ratio = n1/n2)
  g <- g + facet_wrap(~FG_row , ncol=1)
  
  
} 

