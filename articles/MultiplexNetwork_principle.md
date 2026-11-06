# Stochastic Block Models for Multiplex networks: principle

## 1. Preliminaries

This vignette illustrates the use of the `estimateMultiplexSBM` function
and the methods accompanying the R6 class `MultiplexSBM_fit`.

The only package required for the analysis is **sbm**.

``` r

library(sbm)
```

## 2. Multiplex network data

Multiplex networks refer to a collection of networks involving the same
sets of nodes, each network corresponding to a given type of
interaction. Each network may also be referred as a layer. Such a
network is studied in Barbillon et al. ([2017](#ref-multiplexSBM)) where
the nodes are researchers and two networks are provided:

- a network representing the relations of advices between researchers,
- a network representing the exchanges between the laboratories the
  researchers belong to.

In Kéfi ([2016](#ref-multiplexKefi)), the nodes are species and each
network corresponds to a type of ecological relation, namely trophic
links, negative non trophic links and positive non trophic links.

## 3. Stochastic Block models for multiplex networks

### 3.1 General formulation of the model

We give here a few details on Stochastic Block Models (SBM) for
multiplex networks.

Assume that one studies $`n`$ nodes and observes $`L`$ types of links,
resulting into $`L`$ matrices $`(Y^\ell_{ij})_{i,j = 1,\dots,n}`$,
$`\ell = 1, \dots, L`$.

> **The multiplex SBM**
>
> The $`n`$ nodes are divided into $`Q`$ blocks with proportions
> $`\pi = (\pi_1, \dots, \pi_Q)`$, this clustering being encoded into
> the latent variables $`(Z_i)_{i=1,\dots,n}`$ such that $`Z_i = q`$ if
> node $`i`$ belongs to block $`q`$. The relations between nodes $`i`$
> and $`j`$ are driven by the blocks they belong to:
> ``` math
> \begin{align*}
>  (Z_i)_{i} \text{ i.i.d.} \qquad & \mathbb{P}(Z_i = q) = \pi_q, \\
>  (Y^1_{ij}, \dots,Y^L_{ij})_{i \neq j} \text{ indep.} \mid (Z_i)_{i} \qquad & (Y^1_{ij}, \dots,Y^L_{ij}) \mid Z_i = q, Z_j = r \sim \mathcal{F}(\alpha_{qr}),
> \end{align*}
> ```
> where $`\mathcal{F}`$ is a $`L`$ dimensional probability distribution
> adapted to the data.

*Examples*

- If all the networks are binary ($`Y^\ell_{ij} \in \{0,1\}`$) then
  $`\mathcal{F}`$ is a $`L`$-variate Bernoulli distribution and
  ``` math
   \alpha_{qr} = (\alpha^{\omega}_{qr}, \omega \in \{0,1\}^L)
  ```
  with
  $`\alpha^{\omega}_{qr} = \mathbb{P}\left((Y^1_{ij}, \dots,Y^L_{ij}) = \omega \mid Z_i = q, Z_j = r\right) \in [0,1]`$.
- If all the networks are weighted with real values
  ($`Y^\ell_{ij} \in \mathbb{R}`$) then $`\mathcal{F}`$ could be a
  $`L`$-dimensional Gaussian distribution and
  ``` math
   \alpha_{qr} = (\mu_{qr},\Sigma_{qr})
  ```
  with $`\mu_{qr} \in \mathbb{R}^L`$, $`\Sigma_{qr}`$ a $`L \times L`$
  covariance matrix.

Note that these two models assume that the $`L`$ interactions are
dependent **conditionally to the clustering**. The integration against
the latent variables $`Z`$ (clustering) implies that the interactions
between two pairs of nodes are also dependent.

A collection of models are implemented in our package `sbm`. These
models can be classified into two groups: either the layers are
dependent conditionally to $`Z`$ or not.

### 3.2 Dependent and independent layers conditionally to $`Z`$

We consider conditional dependence in a few numbers of models:

- for $`L=2`$ binary networks (Bernoulli);
- for any number $`L`$ of layers with Gaussian multivariate
  distributions but restricted to $`\Sigma_{qr} = \Sigma`$ (same
  covariance in any blocks).

In the binary case, the number of parameters of dependent layers
explodes with the number of layers $`L`$ and the results become
difficult to understand. As a consequence, we also implemented inference
methods for independent layers conditionally to the clusterings for any
number of layers.

> **The multiplex SBM with independent layers**
>
> The clustering $`(Z_i)_{i}`$ is the same as above but, conditionally
> to the clusters, the layers are independent: for any $`\ell,i,j`$,
> ``` math
>  Y^\ell_{ij} \mid Z_i = q, Z_j = r \sim_{ind} \mathcal{F}_\ell(\alpha^\ell_{qr}),
> ```
> where any $`\mathcal{F}_\ell`$ can be Bernoulli, Poisson, Gaussian or
> Zero Inflated Gaussian.

### 3.3 Bipartite multiplex networks

All the previous models can be extended to bipartite networks (i.e. when
the nodes in row are different from the nodes in column). In that case,
two clusterings are introduced.

> **The bipartite multiplex SBM**
>
> The $`n_1`$ row nodes are divided into $`Q_1`$ blocks with proportions
> $`\pi^1`$ and the $`n_2`$ column nodes into $`Q_2`$ blocks with
> proportions $`\pi^2`$, encoded into the row and column clusterings
> $`(Z^1_i)_{i=1,\dots,n_1}`$ and $`(Z^2_j)_{j=1,\dots,n_2}`$:
> ``` math
> \begin{align*}
>  (Z^1_i)_{i} \text{ i.i.d.} \qquad & \mathbb{P}(Z^1_i = q) = \pi^1_q, \\
>  (Z^2_j)_{j} \text{ i.i.d.}, \text{ indep. of } (Z^1_i)_{i} \qquad & \mathbb{P}(Z^2_j = r) = \pi^2_r, \\
>  (Y^1_{ij}, \dots,Y^L_{ij})_{i,j} \text{ indep.} \mid (Z^1_i)_{i}, (Z^2_j)_{j} \qquad & (Y^1_{ij}, \dots,Y^L_{ij}) \mid Z^1_i = q, Z^2_j = r \sim \mathcal{F}(\alpha_{qr}),
> \end{align*}
> ```
> with the same distributions $`\mathcal{F}`$ as above (or independent
> layers
> $`Y^\ell_{ij} \mid Z^1_i = q, Z^2_j = r \sim \mathcal{F}_\ell(\alpha^\ell_{qr})`$).

### 3.4 Inference

The maximization of the likelihood is obtained via a Variational version
of the Expectation-Maximization algorithm. The number of blocks is
chosen via a penalized likelihood criterion (ICL). See Barbillon et al.
([2017](#ref-multiplexSBM)) for details.

Note that, in the case of dependent networks, the inference is performed
via the `blockmodels` package, while the `GREMLINS` package is used in
the case of independent networks. As a consequence, the `estimOptions`
arguments are slightly different.

## 4. Data simulation

The function `sampleMultiplexSBM` supplies a method to simulate
multiplex networks. The argument `dependent` is set to `FALSE` by
default (assuming independent layers conditionally to the clusterings).

We simulate hereafter $`2`$ layers. The corresponding matrices can be
plotted using the function `plotMyMultiplexMatrix`.

``` r

nbLayers <- 2
```

### 4.1 Two independent directed layers: one Bernoulli, one Poisson, two clusters

``` r

Nnodes <- 40
blockProp <- c(.4,.6)
set.seed(1)
connectParam <- list(list(mean=matrix(rbeta(4,.5,.5),2,2)),list(mean=matrix(rexp(4,.5),2,2)))
model <- c("bernoulli","poisson")
type <- "directed"
mySampleMultiplexSBM_PB <-
   sampleMultiplexSBM(
   nbNodes = Nnodes,
    blockProp = blockProp,
   nbLayers = nbLayers,
   connectParam = connectParam,
   model=model,
   dimLabels = c('Individuals'), # generic name of the nodes
   type=type,
   seed = 1)
#> [1] "use of sampleMultipartite"
listSBM_PB <- mySampleMultiplexSBM_PB$listSBM
names(listSBM_PB) <- c("Bernoulli","Poisson")
plotMyMultiplexMatrix(listSBM_PB)
```

![](MultiplexNetwork_principle_files/figure-html/Independant%20Poisson%20bernoulli-1.png)

The resulting object is a list of SBM objects (of the SBM class). To
build such an object from an observed matrix use the function
`defineSBM`.

``` r

listSBM_PB
#> $Bernoulli
#> Fit of a Simple Stochastic Block Model -- bernoulli variant
#> =====================================================================
#> Dimension = ( 40 ) - ( 1 ) blocks.  No  edge covariate(s). No  node covariate(s). 
#> =====================================================================
#> * Useful fields 
#>   $nbNodes, $modelName, $dimLabels, $nbBlocks, $nbCovariates, $nbDyads
#>   $blockProp, $connectParam, $covarParam, $covarList, $covarEffect 
#>   $nbNodesCovariates, $nodesCovariates, $nodesCovariatesParam 
#>   $expectation, $indMemberships, $memberships 
#> * R6 and S3 methods 
#>   $rNetwork, $rMemberships, $rEdges, plot(obj), print(obj), coef(obj) 
#> * Additional fields
#>   $probMemberships, $loglik, $ICL, $storedModels, 
#> * Additional methods 
#>   predict, fitted, $setModel, $reorder 
#> 
#> $Poisson
#> Fit of a Simple Stochastic Block Model -- poisson variant
#> =====================================================================
#> Dimension = ( 40 ) - ( 1 ) blocks.  No  edge covariate(s). No  node covariate(s). 
#> =====================================================================
#> * Useful fields 
#>   $nbNodes, $modelName, $dimLabels, $nbBlocks, $nbCovariates, $nbDyads
#>   $blockProp, $connectParam, $covarParam, $covarList, $covarEffect 
#>   $nbNodesCovariates, $nodesCovariates, $nodesCovariatesParam 
#>   $expectation, $indMemberships, $memberships 
#> * R6 and S3 methods 
#>   $rNetwork, $rMemberships, $rEdges, plot(obj), print(obj), coef(obj) 
#> * Additional fields
#>   $probMemberships, $loglik, $ICL, $storedModels, 
#> * Additional methods 
#>   predict, fitted, $setModel, $reorder
```

### 4.2 Two dependent bipartite Gaussian layers: three row clusters, two column clusters

``` r

blockProp <- list(c(0.3, 0.3, 0.4), c(0.5, 0.5))
Q <- sapply(blockProp, function(p) length(p))
nbNodes <- c(80, 30)
connectParam <- list()
connectParam$mu <- vector("list", nbLayers)
connectParam$mu[[1]] <- matrix(0.1, Q[1], Q[2]) + matrix(c(1, 1, 1, 0,
    1, 0), Q[1], Q[2])
connectParam$mu[[2]] <- matrix(-2, Q[1], Q[2]) + matrix(c(1, 3, 2, 1, 2,
    3), Q[1], Q[2])
connectParam$Sigma <- matrix(c(2, 1, 0.1, 4), nbLayers, nbLayers)
model <- rep("gaussian", 2)
mySampleMultiplexSBM_GG <-
  sampleMultiplexSBM(
     nbNodes = nbNodes,
     blockProp = blockProp,
     nbLayers = nbLayers,
     connectParam = connectParam,
     model = model,
     type = "bipartite",
     dependent = TRUE,
     dimLabels = c('row', 'col'),
     seed = 1)
listSBM_GG <- mySampleMultiplexSBM_GG$listSBM
plotMyMultiplexMatrix(listSBM_GG)
```

![](MultiplexNetwork_principle_files/figure-html/Dependent%20Gaussian%20bipartite-1.png)

### 4.3 Two dependent directed Bernoulli layers: two clusters

``` r

## MultiplexSBM Bernoulli with dependence
Q <- 2
set.seed(94)
P00 <- matrix(runif(Q * Q), Q, Q)
P10 <- matrix(runif(Q * Q), Q, Q)
P01 <- matrix(runif(Q * Q), Q, Q)
P11 <- matrix(runif(Q * Q), Q, Q)
SumP <- P00 + P10 + P01 + P11
P00 <- P00/SumP
P01 <- P01/SumP
P10 <- P10/SumP
P11 <- P11/SumP
connectParam <- list()
connectParam$prob00 <- P00
connectParam$prob01 <- P01
connectParam$prob10 <- P10
connectParam$prob11 <- P11
model <- rep("bernoulli", 2)
type <- "directed"
nbLayers <- 2
Nnodes <- 40
blockProp <- c(0.6, 0.4)
mySampleMultiplexSBM <-
   sampleMultiplexSBM(
     nbNodes = Nnodes,
     blockProp = blockProp,
     nbLayers = nbLayers,
     connectParam = connectParam,
     model = model,
     type = type,
     dependent = TRUE,
     seed = 1)
listSBM_BB <- mySampleMultiplexSBM$listSBM
plotMyMultiplexMatrix(listSBM_BB)
```

![](MultiplexNetwork_principle_files/figure-html/Dependent%20Bernoulli-1.png)

## 5. Inference

We are now able to perform inference on the multiplex network (search of
the “best” number of clusters in terms of ICL). To do so, we have to
chose the model, i.e. the distribution used in each matrix and the
dependence or independence between matrices conditionally to the
clusters.

### 5.1 Two independent directed layers: one Bernoulli, one Poisson

``` r

res_PB <- estimateMultiplexSBM(listSBM_PB)
```

``` r

res_PB$storedModels
#>   indexModel nbParams nbBlocks       ICL    loglik
#> 1          1        9        2 -1926.202 -1894.847
#> 2          2        2        1 -2362.585 -2355.207
```

One can now plot the reorganized matrices and the predicted values.

``` r

plot(res_PB)
```

![](MultiplexNetwork_principle_files/figure-html/plot%20PB-1.png)

``` r

plot(res_PB, type = 'expected')
```

![](MultiplexNetwork_principle_files/figure-html/plot%20PB%20expected-1.png)

or the graphon, i.e. the estimated connectivities alone, each block
being represented by an interval whose length is its estimated
proportion. The two layers following different distributions, the
connectivities of each layer are divided by their maximum to be
represented on the same color scale.

``` r

plot(res_PB, type = 'graphon')
```

![](MultiplexNetwork_principle_files/figure-html/plot%20PB%20graphon-1.png)

One can also compare the estimated clusters to the simulated clusters:

``` r

table(mySampleMultiplexSBM_PB$memberships$Individuals, res_PB$memberships$Individuals)
#>    
#>      1  2
#>   1  0 15
#>   2 25  0
```

We recover the clusters perfectly, up to the label switching of the
blocks.

### 5.2 Two dependent bipartite Gaussian layers

``` r

res_GG <- estimateMultiplexSBM(listSBM_GG, dependent = TRUE, estimOptions = list(plot = FALSE, verbosity = 0 ))
```

``` r

res_GG$storedModels
#>   indexModel nbParams rowBlocks colBlocks nbBlocks       ICL    loglik
#> 5          5       18         3         2        5 -9467.992 -9394.098
#> 6          6       23         4         2        6 -9487.963 -9394.926
#> 7          7       28         5         2        7 -9510.706 -9398.524
#> 8          8       33         6         2        8 -9534.096 -9402.771
#> 4          4       13         2         2        4 -9542.572 -9487.822
#> 3          3        8         2         1        3 -9584.782 -9548.686
#> 2          2        5         1         1        2 -9667.705 -9642.276
#> 1          0       NA         0         0        0        NA        NA
```

``` r

plot(res_GG)
```

![](MultiplexNetwork_principle_files/figure-html/plot%20GG-1.png)

``` r

plot(res_GG, type = 'graphon', plotOptions = list(legend = TRUE))
```

![](MultiplexNetwork_principle_files/figure-html/plot%20GG%20graphon-1.png)

### 5.3 Two dependent directed Bernoulli layers

``` r

res_BB <- estimateMultiplexSBM(listSBM_BB,dependent =  TRUE, estimOptions = list(plot = FALSE, verbosity = 0 ))
```

``` r

res_BB$storedModels
#>   indexModel nbParams nbBlocks       ICL    loglik
#> 2          2       13        2 -1836.174 -1786.056
#> 3          3       29        3 -1890.421 -1778.117
#> 4          4       51        4 -1963.078 -1764.451
#> 1          1        3        1 -2061.409 -2049.340
```

## References

Barbillon, Pierre, Sophie Donnet, Emmanuel Lazega, and Avner Bar-Hen.
2017. “Stochastic block models for multiplex networks: an application to
a multilevel network of researchers.” *Journal of the Royal Statistical
Society Series A* 180 (1): 295–314.
<https://ideas.repec.org/a/bla/jorssa/v180y2017i1p295-314.html>.

Kéfi, Vincent AND Wieters, Sonia AND Miele. 2016. “How Structured Is the
Entangled Bank? The Surprisingly Simple Organization of Multiplex
Ecological Networks Leads to Increased Persistence and Resilience.”
*PLOS Biology* 14 (8): 1–21.
<https://doi.org/10.1371/journal.pbio.1002527>.
