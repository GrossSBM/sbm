## Precompute the fit of the vignette vignettes/Multipartite_EcologicalNetwork.qmd
## and save it in vignettes/Multipartite_EcologicalNetwork_fits.rda, which the
## vignette loads instead of running the fit.
##
## To be run from the root of the package, after any change of the data or
## of the estimation code:
##   Rscript dev/precompute_Multipartite_EcologicalNetwork.R

fitsFile <- "vignettes/Multipartite_EcologicalNetwork_fits.rda"
if (file.exists(fitsFile)) file.remove(fitsFile)

## Knit the vignette: without the rda file, the fit is run
env <- new.env()
knitr::opts_chunk$set(fig.path = file.path(tempdir(), "figure/"))
knitr::knit("vignettes/Multipartite_EcologicalNetwork.qmd", output = tempfile(fileext = ".md"),
            envir = env, quiet = TRUE)

save(list = "myMSBM", envir = env, file = fitsFile)
