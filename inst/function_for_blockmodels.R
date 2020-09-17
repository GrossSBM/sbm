

################################ Plot the BlockModel Network

plotNetBM = function(BMobject,Q){
  membership_name <-  BMobject$membership_name
  a <- extractParamBM(BMobject,Q)$alpha
  b <- extractParamBM(BMobject,Q)$pi
  c <- extractParamBM(BMobject,Q)$piRow
  d <- extractParamBM(BMobject,Q)$piCol

  if (membership_name == 'SBM') {
    G <- graph_from_adjacency_matrix(a, mode = c("undirected"), weighted = TRUE, diag = TRUE)
    plot.igraph(G,vertex.size = b * 100,edge.width = sqrt(abs(E(G)$weight)),vertex.color = 1:Q, layout = layout_nicely)
    # vertex.label=""
  }

  if (membership_name == 'SBM_sym') {
    G <- graph_from_adjacency_matrix(a, mode = c("directed"), weighted = TRUE, diag = TRUE)
    plot.igraph(G,vertex.size = b * 100,edge.width = sqrt(abs(E(G)$weight)),vertex.color = 1:Q, layout = layout_nicely)
    # vertex.label=""
  }

  if (membership_name == 'LBM') {
    G <- graph_from_incidence_matrix(a, weighted = TRUE)
    plot(G,vertex.size = c(c * 100, d * 100), vertex.shape = c("circle", "square")[V(G)$type + 1],
              edge.width = abs(E(G)$weight * 2),vertex.color = 1:Q, layout = layout.bipartite)
  }
  #return(g)
}

#vertex.size=c*100, vertex.size2=d*100,

