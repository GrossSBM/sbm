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
