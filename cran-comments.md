
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
	- Ubuntu Linux 20.04.1 LTS, R-release, GCC

* tested remotely with github-action
  - Linux Ubuntu 22.04, R-release 
  - Linux Ubuntu 22.04, R-oldrel 
  - Linux Ubuntu 22.04, R-devel 
  - Windows Server 2022, R-release, 64 bit
  - macOS Big Sur 11, R-release

- tested remotely with win-builder (R stable and R under development)
  failed on R old release (purrr not available)
  
## R CMD check results (local)

── R CMD check results ────────────────────────────────────────── sbm 0.4.5 ────
Duration: 2m 55s

0 errors ✓ | 0 warnings ✓ | 0 notes ✓

R CMD check succeeded

## R CMD check results (CRAN setup)

For a couple of platforms, a note due to the presence of large components (armadillo library)
