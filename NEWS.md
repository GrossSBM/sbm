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
