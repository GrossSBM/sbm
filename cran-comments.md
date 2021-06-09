
## Version 0.4.3

Minor release that solves a reverse dependency problem i the upcoming version of ggplot2 (v3.3.4)

## Tested environments

* Ubuntu 20.04, R-release GCC (R-hub builder)
* Fedora Linux, R-devel, clang, gfortran (R-hub builder)
* Debian Linux, R-devel, clang (R-hub builder)
* Oracle Solaris 10, x86, 32 bit, R-release  (R-hub builder)
* macOS 10.13.6 High Sierra, R-release, CRAN's setup (R-hub builder)
* Windows Server 2008 R2 SP1, R-release, 32/64 bit (R-hub builder)
* local R installation, R 4.1.0, Ubuntu 20.04
* win-builder (R 4.0.5)
* win-builder (R 4.1.0)
* win-builder (R Under development)
* Windows latest (github-action)
* macOS Catalina 10.15, R-release (github action)
* Linux Ubuntu 20.04, R-release, R-devel (github-action)

## R CMD check results (local)

── R CMD check results ────────────────────────────────────────── sbm 0.4.3 ────
Duration: 2m 15.5s

0 errors ✓ | 0 warnings ✓ | 0 notes ✓

R CMD check succeeded

## R CMD check results (CRAN setup)

For a couple of platforms, a note due to the presence of large components (armadillo library)
