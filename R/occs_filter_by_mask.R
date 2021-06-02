#' Function to filter ocurrence points of Species Temporal Data object (STD object)  by mask.
#' @description This function filters ocurrence points of Species Temporal Data object (STD object)  by mask.
#' @param this.species, Species Temporal Data object see \code{\link[hsi]{sp_temporal_data}}.
#' @param mask A binary raster to filter species occurrence data.
#' @return Returns un  sp.temporal.modeling object (list) con las coordenadas de los puntos de presencia, los anios las observaciones, el path de las capas temporales.
#' @export

occs_filter_by_mask <- function(this.species,mask){
  stopifnot(inherits(this.species, "sp.temporal.modeling"))
  spInMask <- raster::extract(mask,this.species$sp_coords)
  in_mask_index <- which(spInMask==1)

  sp.temp.data.mask <- list(sp_coords = this.species$occs_data[in_mask_index,
                                                               this.species$lon_lat_vars],
                            sp_occs_year = this.species$sp_occs_year[in_mask_index],
                            occs_data = this.species$occs_data[in_mask_index,],
                            lon_lat_vars = this.species$lon_lat_vars ,
                            layers_path_by_year = this.species$layers_path_by_year)
  class(sp.temp.data.mask) <- c("list", "sp.temporal.modeling")

  return( sp.temp.data.mask )
}
