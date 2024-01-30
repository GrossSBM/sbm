#' fungus-tree interaction network
#'
#' This data set provides information about $154$ fungi sampled on $51$ tree species.
#'
#' @format A list with the following entries:
#' * fungi_list list of the fungus species names
#' * tree_list list of the tree species names
#' * fungus_tree binary fungus-tree interactions
#' * tree_tree weighted tree-tree interactions (number of common fungal species two tree species host)
#' * covar_tree covariates associated to pairs of trees (namely genetic, taxonomic and geographic distances)
#'
#' @source Vacher, Corinne, Dominique Piou, and Marie-Laure Desprez-Loustau. "Architecture of an antagonistic tree/fungus network: the asymmetric influence of past evolutionary history." PloS one 3.3 (2008): e1740.
"fungusTreeNetwork"
