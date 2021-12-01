#' Temporal projection of the best model.
#' The function projects the niche model for the d
#' @param this_species a "sp.temp.best.model" object see \code{\link[hsi]{find_best_model}}.
#' @param save_dir Folder where the projected model results will be saved.
#' @param sp_mask A mask to crop the projected area of distribution.
#' @param sp_name The name of the modeled species.
#' @param crs_model The crs of the projection layers.
#' @param plot3d Logical. If the models  have 3 varibles an rgl plot will be shown
#' @return  The output is contained in main directory called "temporal_modeling", it is arranged in subdirectories corresponding to each of the projected years. Each subdirectory has:
#' - Ellipsoid model plot. Only If the number of dimensions to model the niche where 2 or 3.
#' - Binary model of distribution.
#' - A continuos suitability raster.
#' Finally a directory called niche_comparations_results, which has barplots of the percent of area/suitabilty gained (or lost) in each year compared to \code{t_0}.
#' @export
temporal_projection <- function(this_species,
                                projection_years,
                                save_dir,sp_mask,sp_name="",
                                crs_model=NULL, plot3d=FALSE){
  stopifnot(inherits(this_species, "sp.temp.best.model"))
  if(is.null(crs_model)) crs_model <- "+proj=lcc +lat_1=17.5 +lat_2=29.5 +lat_0=12 +lon_0=-102 +x_0=2500000 +y_0=0 +datum=WGS84 +units=m +no_defs"
  #extent_capas <- raster::extent(sp_mask)
  projection_yearsA <- names(this_species$layers_path_by_year)
  projection_yearsU <- as.character(projection_years)
  projection_years <- projection_yearsA[which(projection_yearsA == projection_years)]
  if(length(projection_years) == 0L)
    stop("No environmental data for ", projection_yearsU)
  projection_years <- c(min(projection_yearsA),projection_years)
  projection_vars <- names(this_species$best_model_metadata$centroid)
  ellip <- this_species$best_model_metadata
  # Datos de presencia de la sp en el ambiente
  occs_env <- this_species$env_data_train[,projection_vars]

  to_save_dir <- file.path(save_dir,paste0("temporal_modeling_",sp_name))


  if(!dir.exists(to_save_dir)) dir.create(to_save_dir)
#  print(to_save_dir)
  metadata_dir <- file.path(to_save_dir,
                            "ellipsoid_metadata")
  if(!dir.exists(metadata_dir)) dir.create(metadata_dir)

  # Save ellipsoid metadata
  ellipsoid_metadatafile <-  file.path(metadata_dir,
                                        paste0(sp_name,
                                               "_ellip_metadata.txt"))

  capture.output(ellip ,file = ellipsoid_metadatafile)

  proyections <- lapply(projection_years,function(x){

    layers <- this_species$layers_path_by_year[[x]]
    layers_index <-  unlist(lapply(projection_vars,function(x)
      stringr::str_which(string =layers,pattern = x)))

    #model_layers <- lapply(layers[layers_index],function(x){
    #  r1 <- raster(x)
    #  r1 <- setExtent(r1, extent_capas, keepres=TRUE)
    #  return(r1)
    #})
    #model_layers <- stack(model_layers)
    model_layers <- raster::stack(layers[layers_index])

    sp_model <- ellipsoidfit(data = model_layers,
                             centroid =ellip$centroid,
                             covar =  ellip$covariance,
                             level = this_species$ellipsoid_level,
                             threshold = 0.000001,size = 3,
                             plot = plot3d)

    if(length(ellip$centroid)==3 && plot3d){
      # Presencias de la sp en el ambiente
      rgl::points3d(occs_env,size=10)

      # Ejes del elipsoide

      rgl::segments3d(x = ellip$axis_coordinates[[1]][,1],
                 y = ellip$axis_coordinates[[1]][,2],
                 z = ellip$axis_coordinates[[1]][,3],
                 lwd=3)


      rgl::segments3d(x = ellip$axis_coordinates[[2]][,1],
                 y = ellip$axis_coordinates[[2]][,2],
                 z = ellip$axis_coordinates[[2]][,3],
                 lwd=3)

      rgl::segments3d(x = ellip$axis_coordinates[[3]][,1],
                 y = ellip$axis_coordinates[[3]][,2],
                 z = ellip$axis_coordinates[[3]][,3],
                 lwd=3)

    }


    raster_model <-sp_model$suitRaster*1
    suit_train <- raster::extract(raster_model,
                          this_species$test_data[,this_species$lon_lat_vars])
    min_train_pres <- min(suit_train,na.rm = T)

    if(class(sp_mask)=="RasterLayer"){
      raster_model <- sp_mask*raster_model
    }
    binary_model <- raster_model > min_train_pres



    year_dir <- file.path(to_save_dir,x)
    if(!dir.exists(year_dir)) dir.create(year_dir)

    print(year_dir)


    binary_model_name <- paste0(year_dir,"/",sp_name,"_binary",".tif")
    niche_model_name <-  paste0(year_dir,"/",sp_name,"_niche",".pdf")
    niche_model_map <- paste0(year_dir,"/",sp_name,"_suitability",".tif")
    if(length(ellip$centroid)==3 && plot3d)
      rgl::rgl.postscript(niche_model_name,fmt = "pdf")
    raster::writeRaster(binary_model, binary_model_name,overwrite=T)
    raster::writeRaster(raster_model,niche_model_map,overwrite=T)

    raster::crs(binary_model) <- crs_model
    binary_model_repro <- raster::projectRaster(binary_model,
                                        crs="+proj=longlat +datum=WGS84 +no_defs")
    binary_model_repro <- binary_model_repro > min_train_pres
    distribution_area_model <- distribution_area(binary_model_repro)

    density_model <- density(raster_model@data@values,from=-1,to=2,na.rm=T,n=1000)

    density_model <- data.frame(x= density_model$x,
                                density = density_model$y,
                                area_binary=distribution_area_model)
    names(density_model) <- c(paste0("suitability_",x),
                              paste0("density_",x),
                              paste0("area_binary_",x))
    return(density_model)
  })

  comparations_dir <- file.path(to_save_dir,"niche_comparations_results")
  if(!dir.exists(comparations_dir)) dir.create(comparations_dir)

  density_models <- do.call("cbind.data.frame",proyections)
  density_models_per_year <- file.path(comparations_dir,"density_models_per_year.csv")
  write.csv(density_models,density_models_per_year,
            row.names = FALSE)

  binary_vars <- grep(pattern = "area",x = names(density_models),value = T)
  binary_vars_to_plot <- as.numeric(unique(density_models[,binary_vars]))
  years_barplot <- stringr::str_extract(string = binary_vars,pattern ="[0-9]+" )


  year_to_search <- min(this_species$coords_env_data_all$ID_YEAR)

  projection_years <- list.files(to_save_dir,
                                 recursive = F,pattern = "[0-9]+",
                                 full.names = F)
  projection_years <- as.numeric(projection_years)
  to_compare_dir_index <- which(projection_years == year_to_search)


  suitability_models_path <- list.files(to_save_dir,recursive = T,
                                        pattern = "suitability",full.names = T)
  suitability_models_path <- suitability_models_path[grep("*.tif$",suitability_models_path)]
  suitability_model_init_year <- raster::raster(suitability_models_path[to_compare_dir_index])
  suitability_models_to_compare <- raster::stack(suitability_models_path[-to_compare_dir_index])
  names(suitability_models_to_compare) <- projection_years[-to_compare_dir_index]
  reclass_suitability <- c(-2,- 0.0000001,1 ,- 0.0000001,0, 0, 0, 2, -1)
  #reclass_suitability <- c(-2,-0.000001,1 ,0,0, NA, 0.000001, 2, -1)

  rclmat <- matrix(reclass_suitability, ncol=3, byrow=TRUE)

  suit_models <- lapply(names(suitability_models_to_compare),FUN = function(y){

    suitability_change <- suitability_model_init_year - suitability_models_to_compare[[y]]
    suitability_lost <- raster::reclassify(suitability_change ,rclmat)
    pix_freq <- table(suitability_lost@data@values)
    pix_freq <- (pix_freq /sum(pix_freq ))*100
    change_pix <-  pix_freq[1]- pix_freq[3]
    suit_change_file <- paste0("suitability_change_",year_to_search,y,".tif")
    suit_change_file <- file.path(comparations_dir,suit_change_file)
    suit_lost_file <- paste0("suitability_lost_reclass",year_to_search,y,".tif")
    suit_lost_plot <- paste0("suitability_lost_reclass",year_to_search,y,".pdf")
    suit_lost_plot <- file.path(comparations_dir,suit_lost_plot)
    suit_lost_file <- file.path(comparations_dir,suit_lost_file)

    raster::writeRaster(suitability_change,suit_change_file,overwrite=T)

    raster::writeRaster(suitability_change,suit_lost_file,overwrite=T)

    suit_change_df <- data.frame(comparation=paste0(year_to_search,"vs.",y),
                                 suit_change=change_pix)

    #  1 (color gris) es mejor que antes
    # -1 (color rojo) es peor que antes
    #  0 (color )

    pdf(suit_lost_plot,width = 8,height = 8)
    raster::plot(suitability_lost,col=grDevices::heat.colors(225),
                 main=paste0(year_to_search,"_",gsub("X","",y)))
    dev.off()
    return(suit_change_df)
  })

  dir_final_barplots <-file.path(comparations_dir,
                                 paste0("final_results_",sp_name))
  if(!dir.exists(dir_final_barplots)) dir.create(dir_final_barplots)
  barplot_binary_path <- file.path(dir_final_barplots,"binary_areas_gain_lost.pdf")
  barplot_suitchange_path <- file.path(dir_final_barplots,"suitability_change.pdf")

  pdf(file = barplot_binary_path,width = 8,height = 8)
  barplot(height = binary_vars_to_plot,names.arg=years_barplot )
  dev.off()
  suit_change <- do.call("rbind.data.frame",suit_models)
  suit_change_per_year_file <- file.path(comparations_dir,
                                         "suit_change_per_year.csv")

  pdf(barplot_suitchange_path,width = 8,height = 8)
  barplot(height =as.vector(suit_change$suit_change),names.arg=years_barplot[-1],
          main=paste(years_barplot[1],"vs", "years"))
  dev.off()
  print(dir_final_barplots)
  write.csv(suit_change, suit_change_per_year_file,row.names = F)
  return(suit_change)

}
