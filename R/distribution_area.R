#' Function to compute the total area of distribution of a binary raster.
#' @param binaryMap A binary raster projected in decimal coordinates (longitude a and latitude).
#' @export
distribution_area <- function(binaryMap){
  area1 <- raster::area(binaryMap)
  area_Map <- area1*binaryMap
  area_tot <- raster::cellStats(area_Map,"sum")
  return(area_tot)
}
