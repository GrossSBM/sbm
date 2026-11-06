# Stochastic Block Models for Multiplex networks: application

## 1. Preliminaries

This vignette illustrates the use of the `estimateMultiplexSBM` function
and the methods accompanying the R6 class `MultiplexSBM_fit` on the
`war` data set. The models are described in the vignette *Stochastic
Block Models for Multiplex networks: principle*.

The packages required for the analysis are **sbm** plus some others for
data manipulation and representation:

``` r

library(sbm)
library(igraph)
library(aricode)
library(knitr)
library(patchwork) # to put figures side by side
```

To display the composition of the clusters, we use the following
function, which gives a table with the names of the nodes of each
cluster in a column:

``` r

clusterComposition <- function(memberships, nodeNames) {
  clusters <- split(nodeNames, factor(memberships, levels = 1:max(memberships)))
  maxSize  <- max(lengths(clusters))
  tab <- sapply(clusters, function(x) c(x, rep("", maxSize - length(x))))
  colnames(tab) <- paste("Cluster", seq_along(clusters))
  as.data.frame(tab)
}
```

## 2. Data set

### 2.1 Description of the data

The `war` data set comes in the `sbm` package:

``` r

data("war")
```

This data set contains a list of two networks (`belligerent` and
`alliance`) where the nodes are countries; an edge in the network
`belligerent` means that the two countries have been at war at least
once between years 1816 to 2007; an edge in network `alliance` means
that the two countries have had a formal alliance between years 1816 and
2012. The network `belligerent` has less nodes since countries which
have not been at war at all are not considered.

These two networks were extracted from <https://correlatesofwar.org/>
(see Sarkees and Wayman ([2010](#ref-sarkees2010resort)) for war data,
and Gibler ([2008](#ref-gibler2008international)) for formal alliance).
Version 4.0 was used for war data and version 4.1 for formal alliance.

### 2.2 Data manipulation

Since they don’t have the same size, we choose to only consider nodes
(countries) which were at war with at least one other country. This
corresponds to the first 83 nodes in the Alliance network.

``` r

A <- as_adjacency_matrix(war$alliance, sparse = FALSE)
A <- A[1:83, 1:83]
B <- as_adjacency_matrix(war$belligerent, sparse = FALSE)
```

We can start with a plot of this multiplex network:

``` r

netA <- defineSBM(A, model = "bernoulli", dimLabels = "country")
netB <- defineSBM(B, model = "bernoulli", dimLabels = "country")
listNet <- list(Alliance = netA, Belligerence = netB)
plotMyMultiplexMatrix(listNet)
```

![](Multiplex_allianceNwar_case_study_files/figure-html/plot%20data-1.png)

## 3. Multiplex SBM with independent layers

We run the estimation of this multiplex model. By setting
`dependent = FALSE`, we declare that we consider the two layers to be
independent conditionally on the latent block variables.

> **The multiplex SBM with independent layers**
>
> Denoting $`Y^1`$ the alliance network and $`Y^2`$ the war network, the
> $`n = 83`$ countries are divided into $`Q`$ blocks with proportions
> $`\pi`$ and, for $`\ell = 1, 2`$,
> ``` math
>  Y^\ell_{ij} \mid Z_i = q, Z_j = r \sim_{ind} \mathcal{B}ern(\alpha^\ell_{qr}).
> ```

``` r

MultiplexFitIndep <- estimateMultiplexSBM(listNet, dependent = FALSE,
    estimOptions = list(verbosity = 0))
```

We can retrieve the clustering:

``` r

clust_country_indep <- MultiplexFitIndep$memberships[[1]]
kable(clusterComposition(clust_country_indep, rownames(A)),
      caption = "Country clusters (independent layers)")
```

| Cluster 1 | Cluster 2 | Cluster 3 | Cluster 4 | Cluster 5 |
|:---|:---|:---|:---|:---|
| United States of America | Belize | Bahamas | United Kingdom | Bavaria |
| Canada | El Salvador | Cuba | Netherlands | German Democratic Republic |
|  | Colombia | Haiti | Belgium | Baden |
|  | Ecuador | Dominican Republic | Luxembourg | Wuerttemburg |
|  | Bolivia | Jamaica | France | Austria-Hungary |
|  | Uruguay | Trinidad and Tobago | Spain | Austria |
|  |  | Barbados | Portugal | Czechoslovakia |
|  |  | Dominica | Germany | Slovakia |
|  |  | Grenada | German Federal Republic | Malta |
|  |  | St. Lucia | Poland | Albania |
|  |  | St. Vincent and the Grenadines | Hungary | Croatia |
|  |  | Antigua & Barbuda | Czech Republic | Yugoslavia |
|  |  | St. Kitts and Nevis | Italy | Bosnia and Herzegovina |
|  |  | Mexico | Greece | Cyprus |
|  |  | Guatemala | Norway | Bulgaria |
|  |  | Honduras | Denmark | Moldova |
|  |  | Nicaragua | Iceland | Romania |
|  |  | Costa Rica |  | Russia |
|  |  | Panama |  | Estonia |
|  |  | Venezuela |  | Latvia |
|  |  | Guyana |  | Lithuania |
|  |  | Suriname |  | Ukraine |
|  |  | Peru |  | Belarus |
|  |  | Brazil |  | Armenia |
|  |  | Paraguay |  | Georgia |
|  |  | Chile |  | Azerbaijan |
|  |  | Argentina |  | Finland |
|  |  |  |  | Sweden |
|  |  |  |  | Cape Verde |
|  |  |  |  | Sao Tome and Principe |
|  |  |  |  | Guinea-Bissau |

Country clusters (independent layers) {.table .caption-top
style="width:100%;"}

And we can plot the reorganized adjacency matrices, or the graphon,
i.e. the estimated connectivities $`\widehat{\alpha}^\ell_{qr}`$ alone:
each block is represented by an interval whose length is its estimated
proportion $`\widehat{\pi}_q`$, and the rectangle $`q \times r`$ of the
layer $`\ell`$ is colored according to $`\widehat{\alpha}^\ell_{qr}`$,
the two layers sharing the same color scale. The graphon is close to the
plot of the expectations (`type = "expected"`), which colors each pair
of countries $`(i,j)`$ with
$`\sum_{q,r} \widehat{\tau}_{iq}\widehat{\tau}_{jr}\widehat{\alpha}^\ell_{qr}`$,
$`\widehat{\tau}_{iq}`$ being the estimated probability that country
$`i`$ belongs to block $`q`$: the two plots coincide when the clustering
is clear-cut, but the graphon only depends on the parameters, not on the
nodes.

``` r

plot(MultiplexFitIndep)
```

![](Multiplex_allianceNwar_case_study_files/figure-html/plot%20indep-1.png)

``` r

plot(MultiplexFitIndep, type = "graphon")
```

![](Multiplex_allianceNwar_case_study_files/figure-html/plot%20indep-2.png)

## 4. Multiplex SBM with dependent layers

Now we assume that the two layers are dependent conditionally to the
latent block variables. We then set `dependent = TRUE`.

> **The multiplex SBM with dependent layers**
>
> Conditionally to the blocks, the pair $`(Y^1_{ij}, Y^2_{ij})`$ follows
> a bivariate Bernoulli distribution: for any
> $`\omega = (\omega_1, \omega_2) \in \{0,1\}^2`$,
> ``` math
>  \mathbb{P}\left((Y^1_{ij}, Y^2_{ij}) = \omega \mid Z_i = q, Z_j = r\right) = \alpha^{\omega}_{qr},
> ```
> with
> $`\alpha^{00}_{qr} + \alpha^{01}_{qr} + \alpha^{10}_{qr} + \alpha^{11}_{qr} = 1`$.
> These probabilities are stored in `connectParam$prob00`,
> `connectParam$prob01`, `connectParam$prob10` and
> `connectParam$prob11`.

``` r

MultiplexFitdep <- estimateMultiplexSBM(listNet, dependent = TRUE,
    estimOptions = list(verbosity = 0))
```

We can retrieve the clustering and compare it to the one obtained in the
independent case.

``` r

clust_country_dep <- MultiplexFitdep$memberships[[1]]
kable(clusterComposition(clust_country_dep, rownames(A)),
      caption = "Country clusters (dependent layers)")
```

| Cluster 1 | Cluster 2 | Cluster 3 | Cluster 4 |
|:---|:---|:---|:---|
| United States of America | Bahamas | United Kingdom | Bavaria |
| Canada | Cuba | Netherlands | German Democratic Republic |
|  | Haiti | Belgium | Baden |
|  | Dominican Republic | Luxembourg | Wuerttemburg |
|  | Jamaica | France | Austria-Hungary |
|  | Trinidad and Tobago | Spain | Austria |
|  | Barbados | Portugal | Czechoslovakia |
|  | Dominica | Germany | Slovakia |
|  | Grenada | German Federal Republic | Malta |
|  | St. Lucia | Poland | Albania |
|  | St. Vincent and the Grenadines | Hungary | Croatia |
|  | Antigua & Barbuda | Czech Republic | Yugoslavia |
|  | St. Kitts and Nevis | Italy | Bosnia and Herzegovina |
|  | Mexico | Greece | Cyprus |
|  | Belize | Norway | Bulgaria |
|  | Guatemala | Denmark | Moldova |
|  | Honduras | Iceland | Romania |
|  | El Salvador |  | Russia |
|  | Nicaragua |  | Estonia |
|  | Costa Rica |  | Latvia |
|  | Panama |  | Lithuania |
|  | Colombia |  | Ukraine |
|  | Venezuela |  | Belarus |
|  | Guyana |  | Armenia |
|  | Suriname |  | Georgia |
|  | Ecuador |  | Azerbaijan |
|  | Peru |  | Finland |
|  | Brazil |  | Sweden |
|  | Bolivia |  | Cape Verde |
|  | Paraguay |  | Sao Tome and Principe |
|  | Chile |  | Guinea-Bissau |
|  | Argentina |  |  |
|  | Uruguay |  |  |

Country clusters (dependent layers) {.table .caption-top}

``` r



plotAlluvial(list(independent = clust_country_indep, dependent = clust_country_dep))
```

![](Multiplex_allianceNwar_case_study_files/figure-html/clust%20dep-1.png)

From independent to dependent, the model gathers blocks 3 and 2, i.e put
together all the countries of Latino and south America.

On top of the clustering comparison, we can compare the ICL criterion to
see which of the dependent or independent models is a best fit:

``` r

MultiplexFitdep$ICL
#> [1] -1470.276
MultiplexFitIndep$ICL
#> [1] -1434.602
```

We can do the same plots for the reorganized matrices and the graphon.
Note that the graphon represents the marginal connectivity of each
layer, $`\widehat{\alpha}^{10}_{qr} + \widehat{\alpha}^{11}_{qr}`$ for
the alliances and
$`\widehat{\alpha}^{01}_{qr} + \widehat{\alpha}^{11}_{qr}`$ for the
wars, and it may be relevant to have a look at the conditional
probabilities.

``` r

plot(MultiplexFitdep)
```

![](Multiplex_allianceNwar_case_study_files/figure-html/plot%20dep-1.png)

``` r

plot(MultiplexFitdep, type = "graphon")
```

![](Multiplex_allianceNwar_case_study_files/figure-html/plot%20dep-2.png)

For instance, we may want to compare, for two countries in given blocks,
the probability of having been at war while being allied at some point
(before or after),
$`\alpha^{11}_{qr} / (\alpha^{11}_{qr} + \alpha^{10}_{qr})`$, to the
marginal probability of having been at war,
$`\alpha^{11}_{qr} + \alpha^{01}_{qr}`$:

``` r

p11 <- MultiplexFitdep$connectParam$prob11
p01 <- MultiplexFitdep$connectParam$prob01
p10 <- MultiplexFitdep$connectParam$prob10
round(p11/(p11 + p10), 2)
#>      [,1] [,2] [,3] [,4]
#> [1,] 0.01 0.14 0.03 0.43
#> [2,] 0.14 0.11 0.09 0.10
#> [3,] 0.03 0.09 0.05 0.03
#> [4,] 0.43 0.10 0.03 0.09
```

These two matrices can be represented as graphons, side by side with
**patchwork**, on the same color scale:

``` r

pCond <- p11/(p11 + p10) # conditional probabilities of being at war while having been or will
# be allied
pMarg <- p11 + p01 # marginal probabilities of being at war
limits <- c(0, max(pCond, pMarg))
optGraphon <- list(limits = limits, legend = TRUE, legend.title = TRUE)
gCond <- plotMyGraphon(pCond, MultiplexFitdep$blockProp, nbNodes = 83, dimLabels = "country",
  plotOptions = c(optGraphon, title = "P(war | alliance)"))
gMarg <- plotMyGraphon(pMarg, MultiplexFitdep$blockProp, nbNodes = 83, dimLabels = "country",
  plotOptions = c(optGraphon, title = "P(war)"))
(gCond | gMarg) + plot_layout(guides = "collect") & 
  ggplot2::labs(fill = "Probability") & ggplot2::theme(legend.position = "bottom")
```

![](Multiplex_allianceNwar_case_study_files/figure-html/conditional%20probabilities%20graphon-1.png)

## References

Gibler, Douglas M. 2008. *International Military Alliances, 1648-2008*.
CQ Press.

Sarkees, Meredith Reid, and Frank Whelon Wayman. 2010. *Resort to War: A
Data Guide to Inter-State, Extra-State, Intra-State, and Non-State Wars,
1816-2007*. Cq Pr.
