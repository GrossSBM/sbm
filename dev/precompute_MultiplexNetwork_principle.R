## Precompute the fits of the vignette vignettes/MultiplexNetwork_principle.qmd
## and save them in vignettes/MultiplexNetwork_principle_fits.rda, which the
## vignette loads instead of running the fits.
##
## To be run from the root of the package, after any change of the simulations
## or of the estimation code:
##   Rscript dev/precompute_MultiplexNetwork_principle.R

fitsFile <- "vignettes/MultiplexNetwork_principle_fits.rda"
if (file.exists(fitsFile)) file.remove(fitsFile)

## Knit the vignette: without the rda file, all the fits are run
env <- new.env()
knitr::opts_chunk$set(fig.path = file.path(tempdir(), "figure/"))
knitr::knit("vignettes/MultiplexNetwork_principle.qmd", output = tempfile(fileext = ".md"),
            envir = env, quiet = TRUE)

save(list = c("res_PB", "res_GG", "res_BB"), envir = env, file = fitsFile)
