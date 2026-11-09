## Submission

This is an update of sbm (0.4.7, September 2024) to version 1.0.0.

* It requires blockmodels (>= 1.2.1), which is now on CRAN. blockmodels 1.2.1 fixes the
  predictions of its Bernoulli models with covariates, which makes the tests of sbm 0.4.7
  fail (as reported in the reverse dependency checks of the blockmodels submission). sbm
  1.0.0 fixes its own predictions accordingly.
* New features: covariates on the nodes, missing values in the network, new graphics
  (all plots are now ggplot2 objects), vignettes rewritten with Quarto.
* Bug fixes, some of which change the results (see NEWS.md).
* New authors and contributors.

## Test environments

* local: Ubuntu 24.04, R 4.6.2
* GitHub Actions: ubuntu-latest (R devel, release, oldrel-1), macos-latest (release),
  windows-latest (release)
* win-builder: R devel and release

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

We checked the 3 reverse dependencies on CRAN (missSBM, normalblockr, shinySbm) with
blockmodels 1.2.1, comparing R CMD check results with the CRAN version of sbm (0.4.7) and
with this version: no new problem (all OK in both cases).
