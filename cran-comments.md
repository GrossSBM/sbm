
## Changes since last CRAN release (major)


* support for multipartite SBM
* support for multiplex SBM
* improvement of the plot methods
* simplification/restructuration of classes
* various bug fixes (prediction of Bipartite and Simple SBM, selectModels for Bipartite SBM, estimOptions not taken into account,  display of field covarArray in the absence of covariate)

## Tested environments

* Windows Server 2008 R2 SP1, R-release, 32/64 bit (R-hub builder)
* Ubuntu 20.04, R-release GCC (R-hub builder)
* Fedora Linux, R-devel, clang, gfortran (R-hub builder)
* macOS 10.13.6 High Sierra, R-release, CRAN's setup (R-hub builder)
* Oracle Solaris 10, x86, 32 bit, R-release  (R-hub builder)
* Windows latest (github-action)
* macOS Catalina 10.15, R-release (github action)
* Linux Ubuntu 16.04, R-release (github-action)
* local R installation, R 4.0.5, Ubuntu 20.04
* win-builder (R version old release)
* win-builder (R version release)

## R CMD check results (local)

── R CMD check results ────────────────────────────────────────── sbm 0.4.1 ────
Duration: 2m 15.5s

0 errors ✓ | 0 warnings ✓ | 0 notes ✓

R CMD check succeeded

## R CMD check results (CRAN setup)

Some note due to possibly mis-spelled words in DESCRIPTION (which are false positive)
  al (18:52, 18:101)
  Barbillon (18:39)
  et (18:49, 18:98)
  Multipartite (16:45)
  
