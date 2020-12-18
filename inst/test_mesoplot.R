nbNodes  <- 90
blockProp <- c(.5, .25, .25) # group proportions
means <- diag(.4, 3) + 0.1  # connectivity matrix: affiliation network
# In Bernoulli SBM, parameters is a list with a
# matrix of means 'mean' which are probabilities of connection
connectParam <- list(mean = means)

## Graph Sampling Bernoulli
mySampler <- sampleSimpleSBM(nbNodes, blockProp, connectParam, model = 'bernoulli')
mySampler$plot('meso')
mySampler$rMemberships() # sample new memberships
mySampler$rAdjacency()   # sample new adjacency matrix
plot(mySampler)


## Graph Samplingn Poisson
nbNodes  <- 90
blockProp <- c(.5, .25, .25) # group proportions
means <- diag(15., 3) + 5    # connectivity matrix: affiliation network
# In Poisson SBM, parameters is a list with
# a matrix of means 'mean' which are a mean integer value taken by edges
connectParam <- list(mean = means)

## Graph Sampling
mySampler <- sampleSimpleSBM(nbNodes, blockProp, list(mean = means), model = "poisson")
mySampler$plot('meso',plotOptions = list(vertex.size = 1.4))

#####################################"
library(igraph)



#---------------------------------------------- BIPARTITE

nbNodes <- c(100, 120)
blockProp <- list(c(.5, .5), c(1/3, 1/3, 1/3)) # group proportions
means <- matrix(runif(6), 2, 3)  # connectivity matrix
# In Bernoulli SBM, parameters is a list with
# a matrix of means 'mean' which are probabilities of connection
connectParam <- list(mean = means)
mySampler <- sampleBipartiteSBM(nbNodes, blockProp, connectParam, model = 'bernoulli')


thetaMean <- mySampler$connectParam$mean
pi  <- mySampler$blockProp
directed <- mySampler$directed
bipartite <- TRUE
nbNodes <- mySampler$nbNodes
nodeLabels <- list(row='Book',col = 'Reader')
model  = mySampler$modelName
plotOptions = list()
plotOptions$vertex.size = 1.5
plotMeso(thetaMean, pi,model,directed,bipartite,nbNodes,nodeLabels,plotOptions)



mySampler$plot(type='meso')

