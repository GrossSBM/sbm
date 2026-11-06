## Precompute the fits of the vignette vignettes/SBM_fungus_tree_network.qmd
## and save them in vignettes/SBM_fungus_tree_network_fits.rda, which the
## vignette loads instead of running the fits.
##
## To be run from the root of the package, after any change of the data or
## of the estimation code:
##   Rscript dev/precompute_SBM_fungus_tree_network.R

fitsFile <- "vignettes/SBM_fungus_tree_network_fits.rda"
if (file.exists(fitsFile)) file.remove(fitsFile)

## Knit the vignette: without the rda file, all the fits are run
env <- new.env()
knitr::opts_chunk$set(fig.path = file.path(tempdir(), "figure/"))
knitr::knit("vignettes/SBM_fungus_tree_network.qmd", output = tempfile(fileext = ".md"),
            envir = env, quiet = TRUE)

## Drop the blockmodels object stored in the fits (several MB per fit): the
## results used in the vignette are already imported in the fit. It is kept
## for mySimpleSBM, whose $storedModels and $setModel() are used in the vignette.
for (name in c("mySimpleSBMPoisson", "mySimpleSBMCov", "myBipartiteSBM")) {
  fit <- get(name, env)
  fit$.__enclos_env__$private$BMobject <- NULL
}

save(list = c("mySimpleSBM", "mySimpleSBMPoisson", "mySimpleSBMCov", "myBipartiteSBM"),
     envir = env, file = fitsFile)
