#' Function create a Species Temporal Data object (STD object).
#' @description This function creates Species Temporal Data object.
#' @param occs A data.frame or a SpatialPointsDataFrame with coordinates of the occurrence records. The data must have a time variable indicating the year of each observation.
#' @param longitude If occs is a data.frame the user must indicate the variable name of logitude data.
#' @param latitude If occs is a data.frame the user must indicate the variable name of latitude data.
#' @param sp_year_var A time variable indicating the year of each observation.
#' @param  layers_by_year_dir El directorio que contiene las carpetas donde estan las capas raster que se usaran para modelar; estas deben estar clasificadas por anio.
#' @param layers_ext Extension de los archivos de las capas raster.
#' @param reclass_year_data Logical. If TRUE, occurrence data whose years are not in the environmental layers will be reclassified to the nearest year in the layers. Other wise they will not be used for modeling.
#' @return Returns un  sp.temporal.modeling object (list) con las coordenadas de los puntos de presencia, los anios las observaciones, el path de las capas temporales.
#' @importFrom rgdal readOGR
#' @importFrom raster raster stack
#' @export

sp_temporal_data <- function(occs=NA,longitude=NULL,latitude=NULL,sp_year_var=NA,
                             layers_by_year_dir=NA,layers_ext="*.tif$",
                             reclass_year_data=TRUE){
  classes_sp_data <- c("SpatialPointsDataFrame","data.frame")
  if(class(occs) %in% classes_sp_data){
    if(class(occs) == "SpatialPointsDataFrame"){
      lon_lat_vars <-colnames(occs@coords)
      occs <- data.frame(occs@coords,occs@data)
    }

    if(class(occs) == "data.frame") {
      lon_lat_vars <- c(longitude,latitude)
      if(!all(lon_lat_vars %in% names(occs)))
        stop("\n Please provide species a valid 'longitude' and 'latitude'")
    }
    if(sp_year_var %in% names(occs)){
      # Species year data
      years_in_occs <- as.numeric(as.character(occs[,sp_year_var]))
      years_in_occs_index <- which(!is.na(years_in_occs))
      years_sps <- unique(years_in_occs[years_in_occs_index])
      # environmental layers by year
      layers_all <- list.files(layers_by_year_dir,recursive = T,
                               pattern = layers_ext,full.names = T)

      layers_year_index <- grep(pattern = "[0-9]+",layers_all)
      # years of envrinmental information
      years_layersb <- stringr::str_extract(string = layers_all[layers_year_index],
                                            pattern = "[0-9]+")
      years_layersb <- as.numeric(years_layersb)

      years_layers <- unique(years_layersb)
      layers_path_by_year <- lapply( years_layers,function(x) {
        l_index <- which( years_layersb %in% x)
        layers_all[l_index]
      })
      names(layers_path_by_year) <- years_layers
      # occs data in environmental data
      years_in_layers_data <- years_sps[which(years_sps %in% years_layers)]
      years_not_in_layers_data <- years_sps[-which(years_sps %in% years_layers)]
      if(reclass_year_data){
        to_reclass <- unlist(lapply(years_not_in_layers_data, function(x){
          distance_year <- abs(x-years_in_layers_data)
          year_to_class_index <- which(distance_year== min(distance_year)[1])
          return(years_in_layers_data[year_to_class_index])
        }))

        for(i in 1:length(years_not_in_layers_data)){
          to_reclass_index <- which(years_in_occs %in% years_not_in_layers_data[i])
          years_in_occs[to_reclass_index] <- to_reclass[i]
        }

        occs$years_in_occs <- years_in_occs
        occs <- occs[years_in_occs_index,]

      }
      else{
        years_in_occs_index <- which(years_in_occs %in% years_in_layers_data)
        occs$years_in_occs <- years_in_occs
        occs <- occs[years_in_occs_index,]
      }

      sp_temp_data <- list(sp_coords = occs[,c(lon_lat_vars)],
                           sp_occs_year = occs$years_in_occs,
                           occs_data = occs,
                           lon_lat_vars =lon_lat_vars ,
                           layers_path_by_year = layers_path_by_year)
      graphics::hist(occs$years_in_occs,
                     xlab = "Year",
                     ylab="Number of occs points",
                     main="Occs by year")

      class(sp_temp_data) <- c("list", "sp.temporal.modeling")

      return(sp_temp_data)
    }
    else
      stop("\n Please provide species a valid 'year_var' (it must be numeric)")
  }
  else
    stop("\n Please provide species occurrence data as a 'data.frame' or as a 'SpatialPointsDataFrame'")

}
