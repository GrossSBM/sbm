## Precompute the fits of the vignette vignettes/SBM_nodes_covariates.qmd
## and save them in vignettes/SBM_nodes_covariates_fits.rda, which the
## vignette loads instead of running the fits.
##
## To be run from the root of the package, after any change of the simulations
## or of the estimation code:
##   Rscript dev/precompute_SBM_nodes_covariates.R

fitsFile <- "vignettes/SBM_nodes_covariates_fits.rda"
if (file.exists(fitsFile)) file.remove(fitsFile)

## Knit the vignette: without the rda file, all the fits are run
env <- new.env()
knitr::opts_chunk$set(fig.path = file.path(tempdir(), "figure/"))
knitr::knit("vignettes/SBM_nodes_covariates.qmd", output = tempfile(fileext = ".md"),
            envir = env, quiet = TRUE)

## Drop the blockmodels object stored in each fit (several MB per fit): the
## results used in the vignette are already imported in the fit, only
## $storedModels, $setModel() and $optimize() would need it.
dropBM <- function(obj) {
  if (inherits(obj, c("SimpleSBM_fit", "BipartiteSBM_fit"))) {
    obj$.__enclos_env__$private$BMobject <- NULL
  } else if (is.list(obj)) {
    lapply(obj, dropBM)
  }
  invisible(NULL)
}
for (name in ls(env)) dropBM(get(name, env))

save(list = c("estimOptions",
              "fitCovar", "fitNoCovar", "fitUseless", "Xuseless",
              "fitNoCovarNoisy", "fitCovarNoisy",
              "x3", "covarSets", "fits",
              "networksBip", "fitBip", "fitBipRowOnly"),
     envir = env, file = fitsFile)
