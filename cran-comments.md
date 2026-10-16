
# Next version (in progress, branch `nodescovariates`)

## TO DO before the CRAN submission (then remove this section)

* wait for `blockmodels` with nodes covariates and NA support to be on CRAN
  (CRAN has 1.1.5, the fork `Polarolouis/blockmodels` is 1.2.1)
* remove the `Remotes: Polarolouis/blockmodels` line in DESCRIPTION
* check that `blockmodels (>= 1.2.1)` in DESCRIPTION matches the version actually published
* run `R CMD check --as-cran`

## Changes

* nodes covariates (Simple and Bipartite SBM) and NA support in the network data
* `reorder()`: block proportions `pi` permuted on the blocks (columns), not on the nodes
* Bipartite fit: nodes covariates on one side only can be given as `NULL`
* Simple fit: `setModel()` to the one-block model no longer keeps the previous nodes covariates parameters
* Bipartite `storedModels`: rows always aligned with the index of the stored model
* validation of the `model` argument in a single helper; `defineSBM()` accepts `ZIgaussian` again
* exploratory scripts moved from `inst/` to `dev/` (excluded from the build)

# sbm 0.4.7

* remove spelling.out file causing NOTE on CRAN
* un-naming `$penalty` for bipartite
* updating arXiv links in DESCRIPTION

* tested locally on Ubuntu Linux 24.04 LTS, R-release, GCC

* tested remotely with R-winbuilder 
  - Windows Old release, unstable, unstable

* tested remotely with github-action
  - Linux Ubuntu 22.04, R-release, R-oldrel , R-devel 
  - Windows Server 2022, R-release, 64 bit
  - macOS Big Sur 11, R-release

## R CMD check results (local)

── R CMD check results ────────────────────────────────────────── sbm 0.4.7 ────
Duration: 2m 52s

0 errors ✓ | 0 warnings ✓ | 0 notes ✓

R CMD check succeeded

## R CMD check results (CRAN setup)

For a couple of platforms, a note due to the presence of large components (armadillo library)
