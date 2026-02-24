#function to download the 19 bioclimatic raster layers from Worldclim database
#' get world clim
#'
#' A function to download the 19 bioclimatic variables from the world clim database
#' @return A raster stack of the 19 bioclimatic world clim variables
get_world_clim <- function(path = "~/Repos/Palicourea/Climatic_niche_overlap/") {

  dir.create(path, recursive = TRUE, showWarnings = FALSE)

  bio <- geodata::worldclim_global(
    var = "bio",
    res = 5,
    path = path
  )

  return(bio)
}
