#' Ecological multipartite interaction network
#'
#' Multipartite network of mutualistic interactions between   plants and pollinators, plants and birds  and plants and ants.
#'
#' @format A list a 3 binary  incidence matrices
#' \describe{
#'   \item{Inc_plant_ant}{Interactions between plants (rows) and ants (cols). Matrix with 141 rows and 30 columns}
#'  \item{Inc_plant_bird}{Interactions between plants (rows) and birds (cols). Matrix with141 rows and 46 columns}
#'   \item{Inc_plant_flovis}{Interactions between plants (rows) and pollinators (cols). Matrix with 141 rows and 173 columns}
#'  ...
#' }
#' @source   Dataset  compiled and conducted at Centro de Investigaciones Costeras La Mancha (CICOLMA), located on the central coast of the Gulf of Mexico, Veracruz, Mexico. \url{https://royalsocietypublishing.org/doi/full/10.1098/rspb.2016.1564} \url{https://github.com/lucaspdmedeiros/multi-network_core_removal/tree/master/data}
"multipartiteEcologicalNetwork"
