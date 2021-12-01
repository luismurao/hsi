#' Extract environmental data by year.
#' @description Extract environmental data by year.
#' @param this_species, Species Temporal Data object see \code{\link[hsi]{sp_temporal_data}}.
#' @param layers_pattern A regular expression to indicate the pattern to search in layers files names.
#' @import purrr
#' @export
extract_by_year <- function(this_species,layers_pattern){
  stopifnot(inherits(this_species, "sp.temporal.modeling"))
  years <- sort(unique(this_species$sp_occs_year))
  capasByRes <- lapply(this_species$layers_path_by_year,
                       function(x){
                         x[grep(layers_pattern,x = x)]
                       })

  years_env_fun <- function(x){

    year_env <- x$years_in_occs[1]
    env_layers <- raster::stack(capasByRes[[paste0(year_env)]])
    env_sp <- raster::extract(x = env_layers,
                      y = x[,this_species$lon_lat_vars],
                      df=TRUE)

    names(env_sp) <- c("ID_YEAR",
                       names(env_sp)[2:length(names(env_sp))])
    coords.x1 <- x[,this_species$lon_lat_vars[1]]
    coords.x2 <- x[,this_species$lon_lat_vars[2]]

    env_sp$ID_YEAR <- year_env
    n_70 <- floor(dim(env_sp)[1]*.7)
    env_sp$trian_test <- NA
    train_index <- sample(1:dim(env_sp)[1],n_70)
    env_sp$trian_test[train_index] <- "Train"
    env_sp$trian_test[-train_index] <- "Test"
    return(cbind(coords.x1,coords.x2,env_sp))
  }

  df_occs_year <- this_species$oocs_data
  years_env <- df_occs_year %>%
    split(.$years_in_occs) %>%
    purrr::map_df(~years_env_fun(.x))
  names(years_env)[1:2] <- this_species$lon_lat_vars

  train_data <- which(years_env$trian_test=="Train")
  sp.temp.data.env <- list(sp_coords = years_env[,this_species$lon_lat_vars],
                           coords_env_data_all = years_env,
                           env_data_train = years_env[train_data,
                                                      4:(dim(years_env)[2]-1)],
                           env_data_test = years_env[-train_data,
                                                     4:(dim(years_env)[2]-1)],
                           test_data =years_env[-train_data,],
                           sp_occs_year = years_env$ID_YEAR,
                           oocs_data = this_species$oocs_data,
                           lon_lat_vars = this_species$lon_lat_vars,
                           layers_path_by_year = capasByRes)
  class(sp.temp.data.env) <- c("list", "sp.temporal.modeling","sp.temporal.env")


  return(sp.temp.data.env)

}
