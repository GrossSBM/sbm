
## Changes in 0.2.1 (minor)

Sorry to publish a new version just after yhe first release, but it is to correct
some random error on CRAN during unitary tests

* various bug fixes, especially for directed network in simpleSBM, covariates for BipartiteLBM
* tests for consistency amended to avoid random errors during testthat
* tests for covariates and directed networks
* tests for S3 methods (coef, fitted, predict)

## Tested environments

* local R installation, R 4.0.2, Ubuntu 20.04
* macOS 10.13.6 High Sierra, R-release, CRAN's setup (R-hub builder)
* Windows Server 2008 R2 SP1, R-release, 32/64 bit (R-hub builder)
* Linux Debian, R-release (R-hub builder)
* macOS Catalina 10.15, R-release (github action)
* macOS Catalina 10.15, R-devel (github action)
* Linux ubuntu 16.04, R-release (github-action)
* Windows latest (github-action)
* win-builder (R version 3.6.3)
* win-builder (R version 4.0.2)

## R CMD check results

0 errors | 0 warnings | 0 note
