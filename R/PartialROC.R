#' PartialROC - A modification of the function of the Barve & Barve (2016). ENMGadgets \url{https://github.com/narayanibarve}
#' Function PartialROC generates the area under the curve values using bootstrap method. PartialROC is a model evaluation tool, used for
#' continuous model outputs as compared to binary model outputs. This method is specially used for model trained using presence only data.
#' For more details refer DOI: 10.1016/j.ecolmodel.2007.11.008 and check ENMGadgets \url{https://github.com/narayanibarve}.
#' @param valData - Occurence validation data. Must have 3 columns SpName, Longitude, Latitude.
#' @param PredictionFile - It should be a raster class object of a continuous model output.
#' @param E - Amount of error admissible along the Y-axis, given the requirements and conditions of the study (by default =.05). Value should range between 0 - 1
#' @param RandomPercent - Occurrence points to be sampled randomly from the test data for bootstrapping.
#' @param NoOfIteration - Number of iteration for bootstrapping
#' @param parallel Logical to specify if the computation will be done in parallel. default=TRUE.
#' @return OutputFile will have 4 columns, IterationNo, AUC_at_specified_value, AUC_AT_Random, AUC_Ratio. The first row will always have 0 th interation
#' which is the actual Area Under the Curve without bootstrapping. And the rest of the rows contains auc ratio for all the bootstrap.
#' @importFrom purrr map_df
#' @importFrom future %<-%
#' @importFrom magrittr %>%
#' @importFrom grDevices dev.off pdf
#' @importFrom graphics "barplot" "par" "text"
#' @importFrom stats density lm
#' @importFrom utils capture.output combn write.csv
#' @import future
#' @importFrom stats na.omit
#' @useDynLib hsi
#' @import dplyr
#' @export


PartialROC <- function(valData, PredictionFile, E = 0.05,
                        RandomPercent, NoOfIteration,parallel=FALSE)
{

  test_data <- valData[,2:3]
  continuos_mod <- PredictionFile
  E_percent <- E*100
  boost_percent <- RandomPercent
  n_iter <- NoOfIteration

  if(continuos_mod@data@min == continuos_mod@data@max){
    stop("\nModel with no variability.\n")
  }

  continuos_mod <-round((continuos_mod/raster::cellStats(continuos_mod,
                                                         max)) * 1000)
  test_value <- raster::extract(continuos_mod,test_data)
  test_value <- na.omit( test_value)
  #test_value <- unique(test_value)
  classpixels <- data.frame(raster::freq(continuos_mod))
  classpixels <- data.frame(stats::na.omit(classpixels))

  classpixels <- classpixels  %>% dplyr::mutate_(value= ~rev(value),
                                                 count= ~rev(count),
                                                 totpixperclass = ~cumsum(count),
                                                 percentpixels= ~ totpixperclass/sum(count)) %>%
    dplyr::arrange(value)


  error_sens <- 1-(E_percent/100)
  models_thresholds <- classpixels[,"value"]
  fractional_area <- classpixels[,"percentpixels"]
  n_data <- length(test_value)
  n_samp <- ceiling((boost_percent/100)*(n_data))

  big_classpixels <- matrix(rep(models_thresholds,each=n_samp),
                            ncol=length(models_thresholds))


  calc_aucDF <- function(big_classpixels,fractional_area,
                         test_value,n_data,n_samp,error_sens){

    rowsID <- sample(x = n_data,
                     size = n_samp,
                     replace=TRUE)
    test_value1 <- test_value[rowsID]
    omssion_matrix <-   big_classpixels >  test_value1
    sensibility <- 1 - colSums(omssion_matrix)/n_samp
    xyTable <- data.frame(fractional_area,sensibility)
    less_ID <- which(xyTable$sensibility<=error_sens)
    xyTable <- xyTable[-less_ID,]

    xyTable <- xyTable[order(xyTable$fractional_area,
                             decreasing = F),]

    auc_pmodel <- trap_roc(xyTable$fractional_area,
                           xyTable$sensibility)

    auc_prand <- trap_roc(xyTable$fractional_area,
                          xyTable$fractional_area)
    auc_ratio <- auc_pmodel/auc_prand

    auc_table <- data.frame(auc_pmodel,
                            auc_prand,
                            auc_ratio =auc_ratio )
    return(auc_table)

  }


  if(parallel){

    future::plan(future::multiprocess)
    roc_env <- new.env()
    n_cores <- future::availableCores()
    niter_big <- floor(n_iter/n_cores)
    n_runs <- rep(niter_big,n_cores)
    sum_n_runs <- sum(n_runs)
    n_runs[1] <- n_runs[1] + (n_iter - sum_n_runs)

    for(i in 1:length(n_runs)){
      x <- as.character(i)
      roc_env[[x]] %<-% {
        x1 <- 1:n_runs[i]
        auc_matrix1 <- x1 %>%
          purrr::map_df(~calc_aucDF(big_classpixels,
                                    fractional_area,
                                    test_value,n_data,n_samp,
                                    error_sens))
      }
    }
    partial_AUC <- as.list(roc_env)
    rm(roc_env)
    partial_AUC <- do.call(rbind.data.frame,partial_AUC)
    rownames(partial_AUC) <- NULL
    future::plan(future::sequential)

  }
  else{

    partial_AUC <- 1:n_iter %>%
      purrr::map_df(~calc_aucDF(big_classpixels,
                                fractional_area,
                                test_value,n_data,n_samp,
                                error_sens))

  }
  partial_AUC <- data.frame(NoOfIteration=1:NoOfIteration,partial_AUC)
  return(partial_AUC)


}

