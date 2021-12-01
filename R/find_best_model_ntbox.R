#' Function to find the best n-dimensional ellipsoid model using Partial Roc as a performance criteria.
#' @param this_species, Species Temporal Environment "sp.temporal.env" object see \code{\link[hsi]{extract_by_year}}.
#' @param cor_threshold Threshold valuefrom which it is considered that the correlation is high see \code{\link[hsi]{correlation_finder}}.
#' @param nbg_points Number of background points used to compute partial ROC test. See \code{\link[ntbox]{sample_envbg}} for more details.
#' @param omr_criteria Omission rate used to select best models. See \code{\link[ntbox]{ellipsoid_selection}} for more details.
#' @param ellipsoid_level The proportion of points to be included inside the ellipsoid see \code{\link[hsi]{ellipsoidfit}}.
#' @param nvars_to_fit Number of variables that will be used to model.
#' @param E  Amount of error admissible for Partial Roc test (by default =.05). Value should range between 0 - 1. see \code{\link[hsi]{PartialROC}}
#' @param RandomPercent Occurrence points to be sampled in randomly for the boostrap of the Partial Roc test \code{\link[hsi]{PartialROC}}.
#' @param NoOfIteration Number of iteration for the bootstrapping of the Partial Roc test \code{\link[hsi]{PartialROC}}.
#' @param parallel Logical argument to run computations in parallel. Default TRUE
#' @param n_cores Number of cores to be used in parallelization. Default 4
#' @return A "sp.temp.best.model" object with metadata of the best model given the performance of the Partial Roc test.
#' @export

find_best_model_ntbox <- function(this_species,cor_threshold=0.9,
                                  nbg_points=50000,
                                  omr_criteria =0.1,
                                  ellipsoid_level=0.975,
                                  nvars_to_fit=c(2,3),
                                  E = 0.05,
                                  RandomPercent = 50,
                                  NoOfIteration=1000,
                                  parallel=TRUE,
                                  n_cores=6){
  #ntbox_funcs <- system.file("helpers/helpers_from_ntbox.R",package = "hsi")
  #source(ntbox_funcs)
  stopifnot(inherits(this_species, "sp.temporal.env"))
  n_nas <- floor(dim(this_species$env_data_train)[1]*0.1)
  env_train <- this_species$env_data_train
  env_test <-this_species$env_data_test
  rm_layers <- unlist(sapply( 1:dim(env_train)[2], function(x){
    if(length(which(is.na(env_train[,x]))) > n_nas) return(x)
  } ))

  if(!is.null(rm_layers)){
    env_train <- stats::na.omit(env_train[,-rm_layers])
    env_test <- stats::na.omit(env_test[,-rm_layers])
  }
  cat("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")
  cat("The total number of occurrence records that will be used for model validation is:",
      nrow(env_test), "\n")
  cat("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n")
  numericIDs <- which(sapply(env_train, is.numeric))
  cor_matrix <- stats::cor(env_train[,numericIDs])

  find_cor   <- correlation_finder(cor_mat = cor_matrix,
                                   threshold = cor_threshold,
                                   verbose = F)
  cor_filter <- find_cor$descriptors
  year_to_search <- min(as.numeric(names(this_species$layers_path_by_year)))

  env_layers <- raster::stack(this_species$layers_path_by_year[[paste0(year_to_search)]])
  bg_samp <- sample_envbg(env_layers,nbg = nbg_points)
  ids_fit <- which(nvars_to_fit>length(cor_filter))
  if(length(ids_fit)>0L)
    nvars_to_fit <- nvars_to_fit[-ids_fit]
  if(length(nvars_to_fit) ==0L)
    nvars_to_fit <- 2:length(cor_filter)

  seleccion_mods <- hsi::ellipsoid_selection(env_train = env_train,
                                             env_test = env_test,
                                             env_vars = cor_filter,
                                             nvarstest = nvars_to_fit,
                                             level = ellipsoid_level,
                                             mve = TRUE,
                                             ncores = n_cores,
                                             comp_each = 100,
                                             env_bg = bg_samp,
                                             parallel = parallel,
                                             omr_criteria = omr_criteria,
                                             proc = TRUE,
                                             proc_iter = NoOfIteration,
                                             rseed = 111)
  best_mod_vars <- stringr::str_split(seleccion_mods$fitted_vars,",")
  env_all <- na.omit(this_species$coords_env_data_all[,best_mod_vars[[1]]])
  best_model_metadata <- cov_center(data = env_all,
                                    level = ellipsoid_level,
                                    mve = TRUE,
                                    vars = best_mod_vars[[1]])


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
                             model_selection_results =seleccion_mods,
                             ellipsoid_level =ellipsoid_level)
  class(sp.temp.best.model) <- c("list", "sp.temporal.modeling","sp.temporal.env","sp.temp.best.model")


  return(sp.temp.best.model)

}
