
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
