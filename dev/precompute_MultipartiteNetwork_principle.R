## Precompute the fits of the vignette vignettes/MultipartiteNetwork_principle.qmd
## and save them in vignettes/MultipartiteNetwork_principle_fits.rda, which the
## vignette loads instead of running the fits.
##
## To be run from the root of the package, after any change of the simulations
## or of the estimation code:
##   Rscript dev/precompute_MultipartiteNetwork_principle.R

fitsFile <- "vignettes/MultipartiteNetwork_principle_fits.rda"
if (file.exists(fitsFile)) file.remove(fitsFile)

## Knit the vignette: without the rda file, all the fits are run
env <- new.env()
knitr::opts_chunk$set(fig.path = file.path(tempdir(), "figure/"))
knitr::knit("vignettes/MultipartiteNetwork_principle.qmd", output = tempfile(fileext = ".md"),
            envir = env, quiet = TRUE)

save(list = c("myMSBM", "myMSBM_fixed", "myMSBM_init", "myMSBM_NA"), envir = env,
     file = fitsFile, compress = "xz")
