#' Function to clean duplicated coordinates data
#' @description Clean duplicated longitude and latitude data by year using threshold distance
#'              which is a distance
#' @param this.species, Species Temporal Data object see \code{\link[hdm]{sp_temporal_data}}.
#' @param threshold A numerc value representig the distance (in grades) between coordinates
#'            to be considered as a duplicate.
#' @export

clean_dup_by_year <- function(this.species,threshold){
  stopifnot(inherits(this.species, "sp.temporal.modeling"))
  df_occs_year <- this.species$occs_data

  clean_by_year <- df_occs_year %>%
    split(.$years_in_occs,drop=T) %>%
    purrr::map_df(~clean_dup(data = .x,
                             longitude = this.species$lon_lat_vars[1],
                             latitude = this.species$lon_lat_vars[2],
                             threshold = threshold))

  sp.temp.data.clean <- list(sp_coords = clean_by_year[,this.species$lon_lat_vars],
                             sp_occs_year = clean_by_year$years_in_occs,
                             oocs_data = clean_by_year,
                             lon_lat_vars = this.species$lon_lat_vars ,
                             layers_path_by_year = this.species$layers_path_by_year)
  class(sp.temp.data.clean) <- c("list", "sp.temporal.modeling")


  return(sp.temp.data.clean)
}
. <- NULL
