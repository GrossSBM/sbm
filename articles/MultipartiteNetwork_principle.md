# Stochastic Block Models for Multipartite networks: principle

## 1. Preliminaries

This vignette presents the multipartite stochastic block model and
illustrates, on simulated data, the use of the functions
`sampleMultipartiteSBM` and `estimateMultipartiteSBM` and of the methods
accompanying the R6 class `MultipartiteSBM_fit`. An application to a
real ecological network is given in the vignette *Multipartite
Stochastic Block Models*
([`vignette("Multipartite_EcologicalNetwork", package = "sbm")`](https://grosssbm.github.io/sbm/articles/Multipartite_EcologicalNetwork.md)).

The only package required for the analysis is **sbm**:

``` r

library(sbm)
```

## 2. Multipartite networks

Multipartite networks consist in the joint observation of several
networks implying some common individuals. The individuals (or entities
represented by nodes) at stake are partitioned into groups defined by
their nature, for instance plants, pollinators and birds in an
ecological network. In what follows, these groups are referred to as
**functional groups** (they are the `dimLabels` of the networks in
**sbm**).

A multipartite network is a collection of networks, each of them
involving one or two functional groups:

- a *simple* network represents the relations inside a functional group
  (adjacency matrix, directed or not),
- a *bipartite* network represents the relations between the individuals
  of two functional groups (incidence matrix).

The networks may be binary or weighted, and different networks may have
different types of weights.

## 3. Stochastic Block Models for multipartite networks

### 3.1 The model

Assume that $`S`$ functional groups are at stake, the functional group
$`s`$ containing $`n_s`$ nodes. We denote by $`\mathcal{E}`$ the set of
pairs of functional groups $`(s, s')`$ for which a network $`Y^{ss'}`$
is observed: $`Y^{ss'}_{ij}`$ is the interaction between node $`i`$ of
functional group $`s`$ and node $`j`$ of functional group $`s'`$. If
$`s \neq s'`$, $`Y^{ss'}`$ is an incidence matrix (bipartite network);
$`Y^{ss}`$ is an adjacency matrix, symmetric if the relation inside the
functional group $`s`$ is not oriented.

> **The multipartite SBM**
>
> The $`n_s`$ nodes of functional group $`s`$ are divided into $`Q_s`$
> blocks with proportions $`\pi^s = (\pi^s_1, \dots, \pi^s_{Q_s})`$,
> this clustering being encoded into the latent variables
> $`(Z^s_i)_{i=1,\dots,n_s}`$ such that $`Z^s_i = q`$ if node $`i`$ of
> functional group $`s`$ belongs to block $`q`$. The clusterings of the
> functional groups are independent, and the interactions are driven by
> the blocks of the nodes:
> ``` math
> \begin{align*}
>  (Z^s_i)_{i, s} \text{ indep.} \qquad & \mathbb{P}(Z^s_i = q) = \pi^s_q, \\
>  (Y^{ss'}_{ij})_{(s,s') \in \mathcal{E}, i, j} \text{ indep.} \mid (Z^s_i)_{i,s} \qquad & Y^{ss'}_{ij} \mid Z^s_i = q, Z^{s'}_j = r \sim \mathcal{F}_{ss'}(\alpha^{ss'}_{qr}),
> \end{align*}
> ```
> where $`\mathcal{F}_{ss'}`$ is Bernoulli, Poisson, Gaussian or Zero
> Inflated Gaussian, and may differ from one network to another.

The blocks of a functional group are thus defined by the patterns of
connection of its nodes **in all the networks it is involved in**. This
model is introduced in Bar-Hen et al. ([2022](#ref-multipartite)).

### 3.2 Inference

The parameters are estimated with a variational version of the
Expectation-Maximization algorithm, which also supplies the clustering
of the nodes. The numbers of blocks $`(Q_s)_{s=1,\dots,S}`$ are chosen
with a penalized likelihood criterion (ICL), by a forward-backward
exploration: starting from an initial clustering, blocks are split or
merged, one functional group at a time, as long as the ICL increases.

By default, the initial clusterings are obtained by fitting a simple or
bipartite SBM on each network separately (with the package
**blockmodels**); all the combinations of these clusterings are used as
starting points (option `initBM = TRUE`). An initial clustering can also
be given by the user (option `classifInit`, see Section 7).

The networks may contain missing values (`NA`): the corresponding dyads
are assumed to be missing at random and are ignored in the inference.

The inference is performed by the package **GREMLINS**.

## 4. Data simulation

We simulate a multipartite network involving $`2`$ functional groups,
`A` and `B`, of respective sizes $`n_A = 60`$ and $`n_B = 50`$, divided
into $`3`$ and $`2`$ blocks.

``` r

nbNodes   <- c(60, 50)
blockProp <- list(c(0.16, 0.40, 0.44), c(0.3, 0.7)) # proportions of the blocks in each functional group
```

We observe $`3`$ networks:

- `A`-`B`: continuous weighted interactions (Gaussian),
- `B`-`B`: binary directed interactions (Bernoulli),
- `A`-`A`: counts, undirected (Poisson).

The architecture of the multipartite network is given by a matrix with
one row per network, giving the indices of the functional groups in row
and in column.

``` r

archiMultipartite <- rbind(c(1, 2), c(2, 2), c(1, 1))
model    <- c("gaussian", "bernoulli", "poisson")
directed <- c(NA, TRUE, FALSE) # NA for a bipartite network
C1 <- list(mean = matrix(c(6.1, 8.9, 6.6, 9.8, 2.6, 1.0), 3, 2),
           var  = matrix(c(1.6, 1.6, 1.8, 1.7, 2.3, 1.5), 3, 2))
C2 <- list(mean = matrix(c(0.7, 1.0, 0.4, 0.6), 2, 2))
m3 <- matrix(c(2.5, 2.6, 2.2, 2.2, 2.7, 3.0, 3.6, 3.5, 3.3), 3, 3)
C3 <- list(mean = 0.5 * (m3 + t(m3))) # symmetric, the network being undirected
connectParam <- list(C1, C2, C3)
```

``` r

mySampleMSBM <- sampleMultipartiteSBM(nbNodes, blockProp, archiMultipartite, connectParam,
                                      model, directed, dimLabels = c("A", "B"), seed = 2)
listSBM <- mySampleMSBM$listSBM
trueMemberships <- mySampleMSBM$memberships
```

`listSBM` is a list of networks, as defined by the function `defineSBM`.
The data can be plotted with:

``` r

plotMyMultipartiteMatrix(listSBM)
```

![](MultipartiteNetwork_principle_files/figure-html/plot%20data-1.png)

## 5. Estimation with selection of the numbers of blocks

The numbers of blocks are chosen between the bounds given in
`nbBlocksRange`, a list named with the functional groups (default: from
$`1`$ to $`10`$ blocks in each functional group).

``` r

myMSBM <- estimateMultipartiteSBM(listSBM, estimOptions = list(verbosity = 0, nbCores = 1))
```

The selected numbers of blocks and the agreement of the clusterings with
the simulated ones (adjusted Rand index):

``` r

myMSBM$nbBlocks
#> A B 
#> 3 2
sapply(c("A", "B"), function(s) aricode::ARI(myMSBM$memberships[[s]], trueMemberships[[s]]))
#> A B 
#> 1 1
```

The estimated parameters:

``` r

myMSBM$blockProp
#> $A
#> [1] 0.1333333 0.3500000 0.5166667
#> 
#> $B
#> [1] 0.4 0.6
myMSBM$connectParam[[2]]
#> $mean
#>        [,1]      [,2]
#> [1,] 0.6625 0.3950000
#> [2,] 1.0000 0.6088889
```

All the models visited during the exploration are stored; `$setModel()`
makes one of them the current model.

``` r

myMSBM$storedModels
#>   indexModel nbParams nbBlocks A nbBlocks B nbBlocks       ICL    loglik
#> 1          1       25          3          2        5 -5727.889 -5635.616
#> 2          2       17          2          2        4 -6247.506 -6184.561
#> 3          3       10          1          2        3 -7611.141 -7573.768
```

``` r

plot(myMSBM, type = "expected")
```

![](MultipartiteNetwork_principle_files/figure-html/plot%20fit-1.png)

## 6. Estimation for given numbers of blocks

If the lower and upper bounds of `nbBlocksRange` are equal for all the
functional groups, the parameters are estimated for these numbers of
blocks, without model selection. The bounds may also be equal for some
functional groups only.

``` r

myMSBM_fixed <- estimateMultipartiteSBM(listSBM,
  estimOptions = list(verbosity = 0, nbCores = 1, nbBlocksRange = list(A = c(3, 3), B = c(2, 2))))
```

``` r

myMSBM_fixed$nbBlocks
#> A B 
#> 3 2
```

## 7. Initial clustering

An initial clustering can be given with the option `classifInit`: a list
with one vector of block labels per functional group, named with the
functional groups. The search of the numbers of blocks then starts from
this clustering (`initBM` is ignored). Here, we start from a clustering
with a single block in `A` and $`3`$ blocks in `B`.

``` r

set.seed(1)
myClassifInit <- list(A = rep(1, 60), B = sample(1:3, 50, replace = TRUE))
myMSBM_init <- estimateMultipartiteSBM(listSBM,
  estimOptions = list(verbosity = 0, nbCores = 1, classifInit = myClassifInit))
```

``` r

myMSBM_init$nbBlocks
#> A B 
#> 3 2
myMSBM_init$ICL
#> [1] -5727.889
```

With fixed numbers of blocks, `classifInit` must have these numbers of
blocks: it is then the starting point of the variational EM algorithm.

## 8. Missing data

We remove $`10\%`$ of the dyads of each network, at random
(symmetrically for the undirected network).

``` r

addNA <- function(net, prop) {
  Y <- net$networkData
  isNA <- matrix(runif(length(Y)) < prop, nrow(Y), ncol(Y))
  simple <- inherits(net, "SimpleSBM_fit")
  if (simple && !net$directed) isNA <- isNA | t(isNA)
  Y[isNA] <- NA
  defineSBM(Y, net$modelName, type = ifelse(simple, "simple", "bipartite"),
            directed = net$directed, dimLabels = net$dimLabels)
}
set.seed(3)
listSBM_NA <- lapply(listSBM, addNA, prop = 0.1)
```

``` r

myMSBM_NA <- estimateMultipartiteSBM(listSBM_NA, estimOptions = list(verbosity = 0, nbCores = 1))
```

``` r

myMSBM_NA$nbBlocks
#> A B 
#> 3 2
sapply(c("A", "B"), function(s) aricode::ARI(myMSBM_NA$memberships[[s]], trueMemberships[[s]]))
#> A B 
#> 1 1
```

The method `predict` gives the expected value of each dyad, which can be
used to predict the missing ones. For the binary network `B`-`B`, the
predicted probabilities of the missing dyads are higher when there is
actually an edge:

``` r

pred <- myMSBM_NA$predict()
missingBB <- is.na(listSBM_NA[[2]]$networkData)
tapply(pred[[2]][missingBB], listSBM[[2]]$networkData[missingBB], mean)
#>         0         1 
#> 0.5315599 0.6906959
```

## 9. Parallel computing

The variational EM runs are parallelized with the package **future**. If
no plan has been chosen, a `multisession` plan with `nbCores` workers
(option of `estimOptions`, $`2`$ by default) is used during the
estimation. Another plan can be set before the call:

``` r

future::plan(future::multisession, workers = 4)
myMSBM <- estimateMultipartiteSBM(listSBM)
future::plan(future::sequential)
```

## References

Bar-Hen, Avner, Pierre Barbillon, and Sophie Donnet. 2022. “Block Models
for Generalized Multipartite Networks: Applications in Ecology and
Ethnobiology.” *Statistical Modelling* 22 (4): 273–96.
<https://doi.org/10.1177/1471082X20963254>.
