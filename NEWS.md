# sbm 1.0.0 [2026-10-16]

## Covariates on the nodes

* `estimateSimpleSBM()` and `estimateBipartiteSBM()` gain an argument `nodesCovariates`: covariates on the nodes, driving the block memberships (the block proportions depend on the covariates of each node through a multinomial logistic regression). New field `$nbNodesCovariates`
* new vignette *Stochastic Block Models with covariates on nodes*

## Missing data

* simple and bipartite SBM handle missing values (`NA`) in the network: the missing dyads are not taken into account in the likelihood, nor in the penalty of the ICL
* requires **blockmodels** (>= 1.2.1)

## Graphics: all the plot methods return ggplot2 objects

* the mesoscopic view (`type = 'meso'`) of simple, bipartite and multipartite SBM is now drawn with **ggplot2** (instead of **igraph**): each block is a node whose area is proportional to its proportion, the edges (and loops) are colored and sized according to the connectivity. For multipartite SBM, the functional groups are distinguished by the colors of the nodes, with a free layout or one row per functional group (`layout = 'groups'`); the coordinates of the blocks are returned in the attribute `blockLayout` and can be modified and given back with the option `layout`
* new plot type `type = 'graphon'` for simple, bipartite and multiplex SBM: the block connectivity drawn on the unit square, blocks sized by the block proportions (one panel per layer for multiplex SBM, the marginal connectivity of each layer in the dependent case)
* new function `plotMyGraphon()` to draw the graphon of an SBM from its parameters
* matrix plots (`type = 'data'`): each block gets a background whose opacity grows with its mean value (options `blockBackground`, `block.color`, `blockAlpha.min`, `blockAlpha.max`)
* expected networks (`type = 'expected'`) drawn with the colors of the graphon (option `graphonColors`)
* matrix plots: the labels of the dimensions are now the axis titles (instead of facet strips), missing values (diagonal of simple SBM) are drawn in white, the names of the nodes are aligned with their row/column (option `names.size`)
* multipartite and multiplex matrix plots share the block background and the color scale of simple and bipartite SBM; for multiplex SBM, the layers are written on top of the panels and the nodes at the bottom and on the right
* `plotAlluvial()`: names of the axes kept, stripes colored by the first clustering, clusters reordered to limit the crossings, invisible return

## Multiplex SBM

* the blocks of dependent multiplex fits are reordered as in the independent case (decreasing connectivity summed over the layers)
* `estimateMultiplexSBM(dependent = TRUE)` accepts a named list of layers

## Other changes

* no ICL plot in `$optimize()` when `verbosity = 0`
* vignettes rewritten with **quarto**: new vignette *Stochastic Block Models: principle and illustration*, the vignette *SBM_with_covariates* replaced by the fungus/tree case study, models recalled in callouts with a common notation; the fits are precomputed (scripts in `dev/`)
* pkgdown site with Bootstrap 5

# sbm 0.4.7 [2024-09-16]

* remove spelling.out file causing NOTE on CRAN
* un-naming `$penalty` for bipartite
* updating arXiv links in DESCRIPTION

# sbm 0.4.6 [2024-01-31]

* fixing use of \itemize in documentation 
* remove warnings due to deprecated function in igraph

# sbm 0.4.5 [2023-01-04]

* included a fix for purrr version 1.0.0

# sbm 0.4.4 [2022-08-24]

* small bug fixes
* comply to new CRAN policy for HTML5 in documentation files

# sbm 0.4.3 [2021-06-09]

* minor changes
   - comply with new faceting scale checks in ggplot2
   - avoid testing on Windows when multi-core is not available
   - save output of multiplex vignette to save time during CRAN check

# sbm 0.4.2

* minor changes (encoding)

# sbm 0.4.1

* advanced plotting functions for multilayer/multiplex SBM
* enhanced support for multiplex SBM
* simplification/restructuration of classes
* fix in prediction of Bipartite and Simple SBM (#7)
* fix in selectModels for Bipartite SBM

# sbm 0.4.0 - major [unreleased]

* first support for multiplex SBM
* improvements of multipartite classes

# sbm 0.3.0 - major [unreleased]

* support for multipartite SBM
* improvement of the plot methods
* change of denomination for dataset fungusTreeNetwwork
* fix bug that made estimOptions not to be taken into account in optimize
* fix bug for displaying field covarArray in the absence of any covariate

# sbm 0.2.2 - minor release

* added fields dimLabels and covarArray for class SBM
* added fields nbConnectParam, penalty and entropy  for simple and bipartite SBM (sampler and fit)
* changing roundProduct function to work on list of covariates

# sbm 0.2.1 - minor release

* various bug fixes, especially for directed network in Simple SBM, covariates for Bipartite LBM
* tests for consistency amended to hopefully avoid random errors on CRAN
* tests for covariates and directed networks
* tests for S3 methods (coef, fitted, predict)

# sbm 0.2.0 - first CRAN release

* Support for Simple and Bipartite SBM, undirected and directed with Bernoulli, Poisson and Gaussian emission law of the edges
* Added a `NEWS.md` file to track changes to the package.
