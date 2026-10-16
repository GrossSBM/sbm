## Precompute the fits of the vignette vignettes/Multiplex_allianceNwar_case_study.qmd
## and save them in vignettes/Multiplex_allianceNwar_case_study_fits.rda, which the
## vignette loads instead of running the fits.
##
## To be run from the root of the package, after any change of the data or
## of the estimation code:
##   Rscript dev/precompute_Multiplex_allianceNwar_case_study.R

fitsFile <- "vignettes/Multiplex_allianceNwar_case_study_fits.rda"
if (file.exists(fitsFile)) file.remove(fitsFile)

## Knit the vignette: without the rda file, all the fits are run
env <- new.env()
knitr::opts_chunk$set(fig.path = file.path(tempdir(), "figure/"))
knitr::knit("vignettes/Multiplex_allianceNwar_case_study.qmd", output = tempfile(fileext = ".md"),
            envir = env, quiet = TRUE)

save(list = c("MultiplexFitIndep", "MultiplexFitdep"), envir = env, file = fitsFile)
