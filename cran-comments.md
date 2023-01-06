
Resubmission of sbm 0.4.5, fixing reference to moved URL (CRAN maintainer request)

>  Found the following (possibly) invalid URLs:
     URL: https://www.correlatesofwar.org/ (moved to
https://correlatesofwar.org/)
       From: man/war.Rd
             inst/doc/Multiplex_allianceNwar_case_study.html
       Status: 301
       Message: Moved Permanently

## sbm version 0.4.5

Includes a fix for purrr version 1.0.0 which causes failure on some CRAN platforms  

* tested locally on Ubuntu Linux 22.04 LTS, R-release, GCC

* tested remotely with R-hub 
  - Windows Server 2022, R-devel, 64 bit
  - Fedora Linux, R-devel, clang, gfortran
	- Ubuntu Linux 20.04.1 LTS, R-release, GCC

* tested remotely with github-action
  - Linux ubuntu 22.04, R-release 
  - Linux ubuntu 22.04, R-oldrel 
  - Linux ubuntu 22.04, R-devel 
  - Windows Server 2022, R-release, 64 bit
  - macOS Big Sur 11, R-release

- tested remotely with win-builder (R version 4.1.2, R version 4.1.3)

## Tested environments

* Ubuntu 20.04, R-release GCC (R-hub builder)
* Fedora Linux, R-devel, clang, gfortran (R-hub builder)
* Debian Linux, R-devel, clang (R-hub builder)
* Oracle Solaris 10, x86, 32 bit, R-release  (R-hub builder)
* macOS 10.13.6 High Sierra, R-release, CRAN's setup (R-hub builder)
* Windows Server 2008 R2 SP1, R-release, 32/64 bit (R-hub builder)
* local installation, R 4.2.1, Ubuntu 22.04
* win-builder (R 4.0.5)
* win-builder (R 4.1.0)
* win-builder (R Under development)
* Windows latest (github-action)
* macOS 11 Big Sur, R-release (github action)
* Linux Ubuntu 20.04, R-release, R-devel (github-action)

## R CMD check results (local)

── R CMD check results ────────────────────────────────────────── sbm 0.4.4 ────
Duration: 2m 21s

0 errors ✓ | 0 warnings ✓ | 0 notes ✓

R CMD check succeeded

## R CMD check results (CRAN setup)

For a couple of platforms, a note due to the presence of large components (armadillo library)
