
## sbm version 0.4.6

- fixing use of \itemize in documentation 
- remove warnings due to deprecated function in igraph

* tested locally on Ubuntu Linux 22.04 LTS, R-release, GCC

* tested remotely with R-hub 
  - Windows Server 2022, R-devel, 64 bit
	- Ubuntu Linux 20.04.1 LTS, R-release, GCC

* tested remotely with R-winbuilder 
  - Windows Old release, unstable, unstable

* tested remotely with github-action
  - Linux Ubuntu 22.04, R-release 
  - Linux Ubuntu 22.04, R-oldrel 
  - Linux Ubuntu 22.04, R-devel 
  - Windows Server 2022, R-release, 64 bit
  - macOS Big Sur 11, R-release

## R CMD check results (local)

── R CMD check results ────────────────────────────────────────── sbm 0.4.6 ────
Duration: 2m 56s

0 errors ✓ | 0 warnings ✓ | 0 notes ✓

R CMD check succeeded

## R CMD check results (CRAN setup)

For a couple of platforms, a note due to the presence of large components (armadillo library)
