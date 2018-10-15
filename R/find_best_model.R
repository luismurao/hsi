#' Function to find the best n-dimensional ellipsoid model using Partial Roc as a performance criteria.
#' @param this_species, Species Temporal Environment "sp.temporal.env" object see \code{\link[hsi]{extract_by_year}}.
#' @param cor_threshold Threshold valuefrom which it is considered that the correlation is high see \code{\link[hsi]{correlation_finder}}.
#' @param ellipsoid_level The proportion of points to be included inside the ellipsoid see \code{\link[hsi]{ellipsoidfit}}.
#' @param nvars_to_fit Number of variables that will be used to model.
#' @param plot3d Logical. If the models  have 3 varibles an rgl plot will be shown
#' @param E  Amount of error admissible for Partial Roc test (by default =.05). Value should range between 0 - 1. see \code{\link[hsi]{PartialROC}}
#' @param RandomPercent Occurrence points to be sampled in randomly for the boostrap of the Partial Roc test \code{\link[hsi]{PartialROC}}.
#' @param NoOfIteration Number of iteration for the bootstrapping of the Partial Roc test \code{\link[hsi]{PartialROC}}.
#' @return A "sp.temp.best.model" object with metadata of the best model given the performance of the Partial Roc test.
#' @export

find_best_model <- function(this_species,cor_threshold=0.9,
                            ellipsoid_level=0.975,nvars_to_fit=3,
                            plot3d=FALSE,
                            E = 0.05,
                            RandomPercent = 50,
                            NoOfIteration=1000){
  stopifnot(inherits(this_species, "sp.temporal.env"))
  n_nas <- floor(dim(this_species$env_data_train)[1]*0.1)
  env_train <- this_species$env_data_train

  rm_layers <- unlist(sapply( 1:dim(env_train)[2], function(x){
    if(length(which(is.na(env_train[,x]))) > n_nas) return(x)
  } ))


  env_train <- stats::na.omit(env_train[,-rm_layers])
  cor_matrix <- stats::cor(env_train)
  find_cor   <- correlation_finder(cor_mat = cor_matrix,
                                   threshold = cor_threshold,
                                   verbose = F)
  cor_filter <- find_cor$descriptors
  combinatoria_vars <- combn(length(cor_filter),nvars_to_fit)

  year_to_search <- min(as.numeric(names(this_species$layers_path_by_year)))
  cat("The total number of models to be tested are: ", dim(combinatoria_vars)[2],"...\n\n")

  env_layers <- raster::stack(this_species$layers_path_by_year[[paste0(year_to_search)]])
  modelos <- lapply(1:dim(combinatoria_vars)[2],function(x){
    cat("Doing model: ", x," of ", dim(combinatoria_vars)[2],"\n")

    # Varaibles filtadas por combinatiria de las mas representativas
    vars_model <- cor_filter[combinatoria_vars[,x]]
    ellip <- cov_center(env_train[,vars_model],
                        level = ellipsoid_level ,vars = vars_model)

    # Datos de presencia de la sp en el ambiente
    occs_env <- this_species$env_data_train[,vars_model]

    # Ajuste del modelo de elipsoide

    sp_model <- ellipsoidfit(data = env_layers[[vars_model]],
                             centroid =ellip$centroid,
                             covar =  ellip$covariance,
                             level = ellipsoid_level,
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

    valData <- this_species$test_data[,c(1,2)]
    valData$sp_name <- "sp"
    valData <- valData[,c(3,1,2)]
    p_roc<- PartialROC(valData = valData,
                        PredictionFile = sp_model$suitRaster,
                        E = E,
                        RandomPercent = RandomPercent,
                        NoOfIteration = NoOfIteration)
    p_roc$model <- paste0(x)

    return(list(model = sp_model$suitRaster,
                pRoc=p_roc[,c("auc_ratio","model")],
                metadata=ellip))

  })
  procs <- lapply(1:length(modelos),function(x) {
    proc <- modelos[[x]][[2]]
  })
  procs <- do.call("rbind.data.frame",procs)
  procs$model <- as.factor(procs$model)

  m1 <- lm(auc_ratio ~ model, data = procs)
  model_means <- sapply(levels(procs$model), function(y){
    model_index <- which(procs$model == y)
    media_model <- mean(procs[model_index,1],na.rm=T)
    return(media_model)
  })

  best_model <-names(model_means)[which(model_means==max(model_means))]

  models_meta_data <- lapply(1:length(modelos), function(x){
    matadata <- modelos[[x]][[3]]
  })

  best_model_metadata <- modelos[[as.numeric(best_model)]][[3]]

  sp.temp.best.model <- list(sp_coords = this_species$sp_coords,
                             coords_env_data_all = this_species$coords_env_data_all,
                             env_data_train = this_species$env_data_train,
                             env_data_test = this_species$env_data_test,
                             test_data = this_species$test_data,
                             sp_occs_year = this_species$sp_occs_year,
                             oocs_data = this_species$oocs_data,
                             lon_lat_vars = this_species$lon_lat_vars,
                             layers_path_by_year = this_species$layers_path_by_year,
                             best_model_metadata= best_model_metadata,
                             ellipsoid_level =ellipsoid_level,
                             pROC_table = procs,
                             models_meta_data=models_meta_data)
  class(sp.temp.best.model) <- c("list", "sp.temporal.modeling","sp.temporal.env","sp.temp.best.model")


  return(sp.temp.best.model)

}
