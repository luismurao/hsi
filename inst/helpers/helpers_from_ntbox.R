#' ellipsoid_selection: Performs variable selection for ellipsoid models
#'
#' @description Performs variable selection for ellipsoid models according to omission rates in the environmental space.
#' @param env_train A data frame with the environmental training data.
#' @param env_test A data frame with the environmental testing data. The default is NULL if given the selection process will show the p-value of a binomial test.
#' @param env_vars A vector with the names of environmental variables to be used in the selection process.
#' @param nvarstest A vector indicating the number of variables to fit the ellipsoids during model selection. It is allowed to test models with a different number of variables (i.e. nvarstest=c(3,6)).
#' @param level Proportion of points to be included in the ellipsoids. This parameter is equivalent to the error (E) proposed by Peterson et al. (2008).
#' @param mve A logical value. If TRUE a minimum volume ellipsoid will be computed using
#' the function \code{\link[MASS]{cov.rob}} of the \pkg{MASS} package. If False the covariance matrix of the input data will be used.
#' @param omr_criteria Omission rate criteria. Value of the omission rate allowed for the selection process. Default NULL see details.
#' @param env_bg Environmental data to compute the approximated prevalence of the model. The data should be a sample of the environmental layers of the calibration area.
#' @param parallel The computations will be run in parallel. Default FALSE
#' @param ncores The number of cores that will be used for the parallel process. By default ntbox will use the total number of available cores less one.
#' @param proc Logical if TRUE a partial roc test will be run.
#' @param proc_iter Numeric. The total number of iterations for the partial ROC bootstrap.
#' @param rseed Logical. Whether or not to set a random seed for partial roc bootstrap. Default TRUE.
#' @param comp_each Number of models to run in each job in the parallel computation. Default 100
#' @return A data.frame with 5 columns: i) "fitted_vars" the names of variables that were fitted; ii) "om_rate" omission rates of the model; iii) "bg_prevalence" approximated prevalence of the model see details section; iv) The rank value of importance in model selection by omission rate; v) The rank value by prevalence after if the value of omr_criteria is passed.
#' @details Model selection occurs in environmental space (E-space). For each variable combination the omission rate (omr) in E-space is computed using the function \code{\link[ntbox]{inEllipsoid}}. The results will be ordered by omr and if the user-specified the environmental background "env_bg" an estimated prevalence will be computed and the results will be ordered also by "bg_prevalence".
#'
#' The number of variables to construct candidate models can be specified by the user in the parameter "nvarstest". Model selection will be run in parallel if the user-specified more than one set of combinations and the total number of models to be tested is greater than 500.
#' If given"omr_criteria" and "bg_prevalence", the results will be shown pondering those models that met the "omr_criteria" by the value of "bg_prevalence".
#' For more details and examples go to \code{\link[ntbox]{ellipsoid_omr}} help.
#' @export
#' @import future
#' @author Luis Osorio-Olvera <luismurao@gmail.com>
#' @references Peterson, A.T. et al. (2008) Rethinking receiver operating characteristic analysis applications in ecological niche modeling. Ecol. Modell., 213, 63–72.
#' @examples
#' \dontrun{
#' # Bioclimatic layers path
#' wcpath <- list.files(system.file("extdata/bios",
#'                                 package = "ntbox"),
#'                     pattern = ".tif$",full.names = TRUE)
#' # Bioclimatic layers
#' wc <- raster::stack(wcpath)
#' # Occurrence data for the giant hummingbird (Patagona gigas)
#' pg <- utils::read.csv(system.file("extdata/p_gigas.csv",
#'                                   package = "ntbox"))
#' # Split occs in train and test
#' pgL <- base::split(pg,pg$type)
#' pg_train <- pgL$train
#' pg_test <- pgL$test
#' # Environmental data for training and testing
#' pg_etrain <- raster::extract(wc,pg_train[,c("longitude",
#'                                             "latitude")],
#'                              df=TRUE)
#' pg_etrain <- pg_etrain[,-1]
#' pg_etest <- raster::extract(wc,pg_test[,c("longitude",
#'                                           "latitude")],
#'                             df=TRUE)
#' pg_etest <- pg_etest[,-1]
#'
#' # Non-correlated variables
#' env_varsL <- correlation_finder(cor(pg_etrain),
#'                                        threshold = 0.8,
#'                                        verbose = F)
#' env_vars <- env_varsL$descriptors
#' # Number of variables to fit ellipsoids (3,5,6 )
#' nvarstest <- c(3,5,6)
#' # Level
#' level <- 0.95
#' # Environmental background to compute the appoximated
#' # prevalence in the prediction
#' env_bg <- raster::sampleRandom(wc,10000)
#'
#' # Selection process
#'
#' e_selct <- ellipsoid_selection(env_train = pg_etrain,
#'                                       env_test = pg_etest,
#'                                       env_vars = env_vars,
#'                                       level = level,
#'                                       nvarstest = nvarstest,
#'                                       env_bg = env_bg,
#'                                       omr_criteria=0.07)
#'
#'# Best ellipsoid model for "omr_criteria" and prevalence
#' bestvarcomb <- stringr::str_split(e_selct$fitted_vars,",")[[1]]
#'
#' # Ellipsoid model projection
#'
#' best_mod <- cov_center(pg_etrain[,bestvarcomb],
#'                               mve = T,
#'                               level = 0.99,
#'                               vars = 1:length(bestvarcomb))
#'
#'
#' # Projection model in geographic space
#'
#' mProj <- ellipsoidfit(wc[[bestvarcomb]],
#'                              centroid = best_mod$centroid,
#'                              covar = best_mod$covariance,
#'                              level = 0.99,size = 3)
#'
#' raster::plot(mProj$suitRaster)
#' points(pg[,c("longitude","latitude")],pch=20,cex=0.5)
#'
# Evaluating the model with partial roc of independent data
#' pg_proc <- pROC(continuous_mod = mProj$suitRaster,
#'                        test_data = pg_test[,c("longitude","latitude")],
#'                        n_iter = 1000,
#'                        E_percent = 5,
#'                        boost_percent = 50,parallel = F)

#' print(pg_proc$pROC_summary)
#' }
suppressMessages(library(furrr))
suppressMessages(library(dplyr))
ellipsoid_selection <- function(env_train,env_test=NULL,env_vars,nvarstest,level=0.95,
                                mve=TRUE,env_bg=NULL,omr_criteria,parallel=F,ncores=NULL,
                                comp_each=100,proc=FALSE,
                                proc_iter=100,rseed=TRUE){
  n_vars <- length(env_vars)
  ntest <- sapply(nvarstest, function(x) choose(n_vars,x))
  nmodels <- sum(ntest)
  cat("-----------------------------------------------------------------------------------------\n")

  cat("\t\t**** Starting model selection process ****\n")
  cat("-----------------------------------------------------------------------------------------\n\n")
  for(i in 1:length(ntest)){
    cat("A total number of",ntest[i] ,"models will be created for combinations",
        "of",n_vars, "variables taken by",nvarstest[i],"\n\n")
  }
  cat("-----------------------------------------------------------------------------------------\n")
  cat("\t **A total number of",nmodels ,"models will be tested **\n\n")
  cat("-----------------------------------------------------------------------------------------\n")



  if(nmodels >100 && parallel){
    max_var <- max(nvarstest)
    cvars <- lapply(nvarstest, function(x) {

      cb <- utils::combn(env_vars,x)
      if(x < max_var){
        nrowNA <-max_var-nrow(cb)
        na_mat <- matrix(nrow = nrowNA,ncol=ncol(cb))
        cb <- rbind(cb,na_mat)
      }
      return(cb)
    })
    big_vars <- do.call(cbind,cvars)

    n_cores <- future::availableCores() -1
    if(ncores>n_cores || is.null(ncores)){
      n_cores <- n_cores
    } else{
      n_cores <- ncores
    }
    niter_big <- floor(nmodels/n_cores)
    if(niter_big>comp_each)
      niter_big <- comp_each
    steps <- seq(1, nmodels, niter_big)
    nsteps <- length(steps)
    if(steps[nsteps]<nmodels){
      kkk <- c(steps,  nmodels + 1)
    } else {
      kkk <- steps
      kkk[nsteps] <- kkk[nsteps] + 1
    }

    long_k <- length(kkk)
    pasos <- 1:(length(kkk) - 1)
    pasosChar <- paste0(pasos)
    globs <- c("env_train",
               "env_test",
               "env_bg")
    #source()
    furrr::furrr_options(globals = c("env_train",
                                     "env_test",
                                     "env_bg",
                                     "rseed","level"),
                         packages = c("Rcpp","ntbox"))
    plan(multisession,workers=n_cores)
    options(future.globals.maxSize= 8500*1024^2)
    model_select <- new.env()
    for (paso in pasosChar) {
      x <- as.numeric(paso)
      #fname <- file.path(dir1,paste0("eselection_",x,".txt"))
      #if(x>n_cores) core <- 1

      cat("Doing calibration from model ",
          kkk[x],"to ",kkk[x + 1] - 1,
          "in process ",x,"\n\n")
      model_select[[paso]] %<-% {

        seq_model <- kkk[x]:(kkk[x + 1] - 1)
        combs_v <- as.matrix(big_vars[,seq_model])

        results_L <- lapply(1:ncol(combs_v),function(x_comb) {
          var_comb <- stats::na.omit(combs_v[,x_comb])
          env_data0 <- stats::na.omit(env_train[,var_comb])
          env_test0 <- stats::na.omit(env_test[,var_comb])
          env_bg0 <-   stats::na.omit(env_bg[,var_comb])
          r1 <- ellipsoid_omr(env_data = env_data0,
                                     env_test = env_test0,
                                     env_bg = env_bg0,
                                     cf_level = level,
                                     proc = proc,
                                     proc_iter,rseed=rseed)
          return(r1)
        })
        results_df <- do.call("rbind.data.frame",results_L)
        cat("Finishing calibration of models ",kkk[x],"to ",kkk[x + 1] - 1,
            "\n\n")
        return(results_df)
      }

    }
    mres <- as.list(model_select)

    cat("Finishing...\n\n")
    cat("-----------------------------------------------------------------------------------------\n")


    rfinal <- do.call("rbind.data.frame", mres )

    future::plan(sequential)
  }
  else{
    cvars <- lapply(nvarstest, function(x) utils::combn(env_vars,x))

    results_L <- lapply(1:length(cvars), function(x) {
      combs_v <- cvars[[x]]
      results_L <- lapply(1:ncol(combs_v),function(x_comb) {
        var_comb <- stats::na.omit(combs_v[,x_comb])
        env_data <- stats::na.omit(env_train[,var_comb])
        env_test <- stats::na.omit(env_test[,var_comb])
        env_bg <-   stats::na.omit(env_bg[,var_comb])
        r1 <- ellipsoid_omr(env_data = env_data,
                                   env_test = env_test,
                                   env_bg = env_bg,
                                   cf_level = level,
                                   proc = proc,
                                   proc_iter,rseed=rseed)
        return(r1)
      })
      results_df <- do.call("rbind.data.frame",results_L)
      return(results_df)
    })
    rfinal <- do.call("rbind.data.frame",results_L)
  }
  bg_omr <- c("bg_prevalence","om_rate_test") %in% names(rfinal)
  bg_omr_in <- all(bg_omr)
  if( bg_omr_in){
    mean_omr <- rowMeans(rfinal[,c("om_rate_train",
                                   "om_rate_test")])
    rfinal$mean_omr_train_test <- mean_omr
    rfinal <- rfinal[order(rfinal$mean_omr_train_tes,
                           rfinal$bg_prevalence,
                           decreasing = F),]

    rfinal <- data.frame(rfinal,rank_by_omr_train_test=1:nrow(rfinal))
    met_criteriaID_train <- which(rfinal$om_rate_train <= omr_criteria)
    met_criteriaID_test <- which(rfinal$om_rate_test <= omr_criteria)
    met_criteriaID_both <- intersect(met_criteriaID_train,
                                     met_criteriaID_test)

    if(length(met_criteriaID_train) > 0L){
      cat("\t",length(met_criteriaID_train),
          "models passed omr_criteria for train data\n")
    }
    if(length(met_criteriaID_test) > 0L){
      cat("\t",length(met_criteriaID_test),
          "models passed omr_criteria for test data\n")

    }
    if(length(met_criteriaID_both) > 0L){
      cat("\t",length(met_criteriaID_both),
          "models passed omr_criteria for train and test data\n")
    }
    else{
      cat("\tNo model passed the omission criteria ranking by mean omission rates\n")
      return(rfinal)
    }
    best_r <- rfinal[met_criteriaID_both,]
    if(proc){
      best_r <- best_r[order(best_r$env_bg_paucratio,
                             decreasing = TRUE),]
    }

    rfinal <- rbind(best_r,
                    rfinal[-met_criteriaID_both,])
    if(proc){
      rfinal <- data.frame(rfinal,
                           rank_omr_aucratio=1:nrow(rfinal))
    }
  }
  else
    rfinal <- rfinal[order(rfinal$om_rate_train,
                           decreasing = F),]
  rownames(rfinal) <- NULL
  return(rfinal)
}


#' ellipsoid_omr
#'
#' @description Compute the omission rate of ellipspoid models
#' @param env_data A data frame with the environmental data.
#' @param env_test A data frame with the environmental testing data. The default is NULL if given the selection process will show the p-value of a binomial test.
#' @param env_bg Environmental data to compute the approximated prevalence of the model. The data should be a sample of the environmental layers of the calibration area.
#' @param cf_level Proportion of points to be included in the ellipsoids. This parameter is equivalent to the error (E) proposed by Peterson et al. (2008).
#' @param mve A logical value. If TRUE a minimum volume ellipsoid will be computed using
#' the function \code{\link[MASS]{cov.rob}} of the \pkg{MASS} package. If False the covariance matrix of the input data will be used.
#' @param proc Logical if TRUE a partial roc test will be run.
#' @param proc_iter Numeric. The total number of iterations for the partial ROC bootstrap.
#' @param rseed Logical. Whether or not to set a random seed for partial roc bootstrap. Default TRUE.
#' @return A data.frame with 5 columns: i) "fitted_vars" the names of variables that were fitted; ii) "om_rate" omission rates of the model; iii) "bg_prevalence" approximated prevalence of the model see details section.
#' @export
#' @examples
#' \dontrun{
#' # Bioclimatic layers path
#' wcpath <- list.files(system.file("extdata/bios",
#'                                 package = "ntbox"),
#'                     pattern = ".tif$",full.names = TRUE)
#' # Bioclimatic layers
#' wc <- raster::stack(wcpath)
#' # Occurrence data for the giant hummingbird (Patagona gigas)
#' pg <- utils::read.csv(system.file("extdata/p_gigas.csv",
#'                                   package = "ntbox"))
#' # Split occs in train and test
#' pgL <- base::split(pg,pg$type)
#' pg_train <- pgL$train
#' pg_test <- pgL$test
#' # Environmental data for training and testing
#' pg_etrain <- raster::extract(wc,pg_train[,c("longitude",
#'                                             "latitude")],
#'                              df=TRUE)
#' pg_etrain <- pg_etrain[,-1]
#' pg_etest <- raster::extract(wc,pg_test[,c("longitude",
#'                                           "latitude")],
#'                             df=TRUE)
#' pg_etest <- pg_etest[,-1]
#'
#' # Non-correlated variables
#' env_varsL <- correlation_finder(cor(pg_etrain),
#'                                        threshold = 0.8,
#'                                        verbose = F)
#' env_vars <- env_varsL$descriptors
#' env_bg <- raster::sampleRandom(wc,10000)
#' ellip_eval <- ellipsoid_omr(env_data=pg_etrain[,c("bio01","bio07","bio12")],
#'                             env_test=pg_etest[,c("bio01","bio07","bio12")],
#'                             env_bg = env_bg[,c("bio01","bio07","bio12")],
#'                             cf_level = 0.97,
#'                             mve=TRUE,proc=TRUE,
#'                             proc_iter=100,rseed=TRUE)
#' print(ellip_eval)
#' }
ellipsoid_omr <- function(env_data,env_test=NULL,env_bg,cf_level,mve=TRUE,proc=FALSE,proc_iter=100,rseed=TRUE){
  emd <- try(cov_center(data = env_data,
                               mve = mve,
                               level = cf_level,
                               vars = 1:ncol(env_data)),
             silent = TRUE)

  message1 <- attr(emd,"class")== "try-error"
  if(length(message1)>0L)
    return()

  in_e <-  inEllipsoid(centroid = emd$centroid,
                       eShape = emd$covariance,
                       env_data = env_data,
                       level = cf_level)

  fails_train_ids <- which(in_e$in_Ellipsoid== 0)

  if(length(fails_train_ids)>0){
    fails_train_ids <- paste0(fails_train_ids,collapse = ",")
  } else {
    fails_train_ids <- NA
  }

  occs_table <- table( in_e$in_Ellipsoid)

  succsID <- which(names(occs_table) %in% "1")
  failsID <- which(names(occs_table) %in% "0")

  occs_succs <-  if(length(succsID)>0L){
    occs_table[[succsID]]
  } else{
    0
  }
  occs_fail <-  if(length(failsID)>0L){
    occs_table[[failsID]]
  } else{
    0
  }

  a_train <-  occs_fail
  omrate_train <- a_train /nrow( in_e)

  d_results <- data.frame(fitted_vars =paste(names(emd$centroid),
                                             collapse =  ","),
                          nvars=length(emd$centroid),
                          om_rate_train= omrate_train,
                          non_pred_train_ids = fails_train_ids)
  if(is.data.frame(env_test) || is.matrix(env_test)){
    in_etest <-  inEllipsoid(centroid = emd$centroid,
                                    eShape = emd$covariance,
                                    env_data = env_test,
                                    level = cf_level)

    fails_test_ids <- which(in_etest$in_Ellipsoid== 0)

    if(length(fails_train_ids)>0){
      fails_test_ids <- paste0(fails_test_ids,collapse = ",")
    } else {
      fails_test_ids <- NA
    }

    suits_val <- exp(-0.5*( in_etest$mh_dist))

    occs_table_test <- table(in_etest$in_Ellipsoid)

    succsID <- which(names(occs_table_test) %in% "1")
    failsID <- which(names(occs_table_test) %in% "0")

    occs_succs_test <-  if(length(succsID)>0L){
      occs_table_test[[succsID]]
    } else{
      0
    }
    occs_fail_test <-  if(length(failsID)>0L){
      occs_table_test[[failsID]]
    } else{
      0
    }
    a_test <-  occs_fail_test
    omrate_test <- a_test /nrow( in_etest)
    d_results <- data.frame(d_results,
                            om_rate_test=omrate_test,
                            non_pred_test_ids=fails_test_ids)
  }

  if(!is.null(env_bg)){

    env_bg <- data.frame(env_bg)
    in_ebg <-  inEllipsoid(centroid = emd$centroid,
                                  eShape = emd$covariance,
                                  env_data = env_bg,
                                  level = cf_level)
    suits_bg <- exp(-0.5*in_ebg$mh_dist)

    bg_table <- table(c(in_ebg$in_Ellipsoid,in_e$in_Ellipsoid))
    succs_bg_ID <- which(names(bg_table) %in% "1")
    fails_bg_ID <- which(names(bg_table) %in% "0")

    bg_succs <-  if(length(succs_bg_ID)>0L){
      bg_table[[succs_bg_ID]]
    }
    else{
      0
    }

    bg_fails <-  if(length(fails_bg_ID)>0L){
      bg_table[[fails_bg_ID]]
    }
    else{
      0
    }
    prevBG <- bg_succs/(bg_fails+bg_succs)
    d_results <-data.frame( d_results,
                            bg_prevalence= prevBG)

    if(exists("in_etest")){
      #bin_table <- table(c(in_ebg$in_Ellipsoid,
      #                     in_etest$in_Ellipsoid))
      #binBG <- bin_table[[2]]/(bin_table[[1]]+bin_table[[2]])
      test_fail <-  occs_fail_test
      test_succs <- occs_succs_test
      p_bin <- 1 - stats::pbinom(test_succs,
                                 size=test_succs+test_fail,
                                 prob = prevBG)
      d_results <-data.frame( d_results,
                              pval_bin=p_bin)
      if(proc){
        proc <- pROC(suits_bg,test_data = suits_val,
                     n_iter = proc_iter,rseed = rseed,E_percent = 5)
        pval_proc <- proc$pROC_summary[3]
        mean_aucratio <- proc$pROC_summary[2]
        mean_auc <- proc$pROC_summary[1]
        d_results <-data.frame( d_results,
                                pval_proc,
                                env_bg_paucratio= mean_aucratio,
                                env_bg_auc = mean_auc)
      }


    }

  }
  return(d_results)
}


#' inEllipsoid: Determine if a point is inside or outside an ellipsoid
#'
#' @description Determine if a point is inside or outside an ellipsoid.
#' @param centroid A numeric vector of centroids for each environmental variable
#' @param eShape Shape matrix of the ellipsoid (can be a covariance matrix or a minimum volume ellipsoid).
#' @param env_data A data frame with the environmental training data.
#' @param level Proportion of points to be included in the ellipsoids. This parameter is equivalent to the error (E) proposed by Peterson et al. (2008).
#' @return A data.frame with 2 columns. The first "in_Ellipsoid" binary response with values 1 (inside the ellipsoid) and zeros (outside the ellipsoid); the second "mh_dist" Mahalanobis distance to centroid.
#' @export
#' @examples
#' \dontrun{
#' # Bioclimatic layers path
#' wcpath <- list.files(system.file("extdata/bios",
#'                                 package = "ntbox"),
#'                     pattern = ".tif$",full.names = TRUE)
#' # Bioclimatic layers
#' wc <- raster::stack(wcpath)
#' # Occurrence data for the giant hummingbird (Patagona gigas)
#' pg <- utils::read.csv(system.file("extdata/p_gigas.csv",
#'                                   package = "ntbox"))
#' # Environmental data
#' pg_env <- raster::extract(wc,pg[,c("longitude",
#'                                       "latitude")],
#'                              df=TRUE)
#' pg_env <- pg_env[,-1]
#'
#' pg_ellip <- cov_center(pg_env,mve=TRUE,
#'                        level=0.95,
#'                        vars = c("bio05",
#'                                 "bio06",
#'                                 "bio12"))
#' # Environmental random data
#' env_rdata <- raster::sampleRandom(wc,1000)
#' inErdata <- inEllipsoid(env_data = env_rdata[,c("bio05",
#'                                                 "bio06",
#'                                                 "bio12")],
#'                         centroid = pg_ellip$centroid,
#'                         eShape=pg_ellip$covariance,
#'                         level = 0.99)
#'
#' }
inEllipsoid <- function(centroid,eShape,env_data,level){

  mh_dist <- stats::mahalanobis(env_data,
                                center = centroid,
                                cov =eShape)
  in_Ellipsoid <- mh_dist <= stats::qchisq(level,
                                           length(centroid))
  in_Ellipsoid <- in_Ellipsoid*1
  in_Ellipsoid_mh <- data.frame(in_Ellipsoid,mh_dist )

  return(in_Ellipsoid_mh)
}


#' Generate random environmental background data
#'
#' @description Generate environmental background data is a function similar
#' to sampleRandom function of the raster package but optimized for Ecological
#' niche modeling.
#' @param envlayers A raster stack or brick.
#' @param nbg Number of points for the background data
#' @param nprop Proportion of environmental data to be sampled. Default NULL
#' @param coordinates Logical. If TRUE cell coordinates will be returned
#' @param cellIDs Logical. If TRUE cell IDs will be returned
#' @param rseed Random seed number. Default NULL
#' @param ncores Number of workers to run the parallel process.
#' @import future
#' @examples
#' \dontrun{
#' wcpath <- list.files(system.file("extdata/bios",
#'                                  package = "ntbox"),
#'                      pattern = ".tif$",
#'                      full.names = TRUE)
#'
#' envlayers <- raster::stack(wcpath)
#' vals <- sample_envbg(envlayers,nbg = 3583)
#' # Using a proportion of data
#' vals <- sample_envbg(envlayers,nprop = 0.20)
#' }
#' @export
sample_envbg <- function(envlayers,nbg,nprop=NULL,coordinates=FALSE,
                         cellIDs=FALSE,rseed=NULL,ncores=4){
  if(class(envlayers) == "RasterStack" ||
     class(envlayers) == "RasterBrick"){
    envlayers <- raster::stack(envlayers)
    l1 <- envlayers[[1]]
    #nona <- raster::Which(!is.na(l1),cells=TRUE)
    nona <- which(!is.na(as.vector(l1)))
    n_nona <- length(nona)
    if(!is.null(nprop)){
      npoints <- ceiling(nprop*n_nona)
    }
    else{
      npoints <- nbg
    }
    if(!is.numeric(rseed))
      set.seed(rseed)
    #cat("Number of points to be sampled:",npoints)
    toSamp <- sample(nona,size = npoints,replace = FALSE)
    canP <- raster::canProcessInMemory(l1,
                                       n=raster::nlayers(envlayers))
    if(canP){
      env_bg <- envlayers[toSamp]
    }
    else {
      n_cores <- future::availableCores() -1
      if(ncores>n_cores || is.null(ncores)){
        n_cores <- n_cores
      } else{
        n_cores <- ncores
      }
      fnames <- sapply(envlayers@layers, function(x) x@file@name)
      fnames <- unique(fnames)
      indexL <- 1:raster::nlayers(envlayers)
      furrr::furrr_options(globals = c("fnames",
                                       "toSamp",
                                       "indexL"))
      plan(multisession,workers=n_cores)
      options(future.globals.maxSize= 8500*1024^2)
      env_bg <- furrr::future_map_dfc(indexL, function(x){
        if(length(fnames) == 1)
          r1 <- raster::raster(fnames,band=x)
        else
          r1 <- raster::raster(fnames[x])
        r2 <- r1[]
        d1 <- data.frame(r2[toSamp])
        names(d1) <- names(r1)
        return(d1)
      },.progress = TRUE)
      future::plan(future::sequential)
    }
    if(coordinates){
      coords <- raster::xyFromCell(l1,toSamp)
      env_bg <- data.frame(coords,env_bg)
    }
    if(cellIDs){
      env_bg <- data.frame(cellID=toSamp,env_bg)
    }

  }
  else
    stop("envlayers should be of class RasterStack or RasterBrick")
  return(env_bg)
}

#' Partial ROC calculation for Niche Models
#'
#' @description pROC applies partial ROC tests to continuous niche models.
#'
#' @param continuous_mod a RasterLayer or a numeric vector of the ecological niche model to be evaluated. If a numeric vector is provided it should contain the values of the predicted suitability.
#' @param test_data A numerical matrix, data.frame, or a numeric vector. If it is data.frame or matrix it should contain coordinates of the occurrences used to test the ecological niche model to be evaluated; columns must be: longitude and latitude. If numeric vector it should contain the values of the predicted suitability.
#' @param E_percent (numeric) value from 0 to 100 that will be used as a threshold (E);
#' default = 5.
#' @param boost_percent (numeric) value from 0 to 100 representing the percent of testing data
#' to be used for performing the bootstrap process for calculating the partial ROC;
#' default = 50.
#' @param n_iter (numeric) number of bootstrap iterations to be performed;
#' default = 1000.
#' @param rseed Logical. Whether or not to set a random seed. Default FALSE.
#' @param parallel Logical to specify if the computation will be done in parallel. default=TRUE.
#' @param ncores Numeric; the number of cores to be used for parallelization.
#' @return A data.frame containing the AUC values and AUC ratios calculated for each iteration.
#' @details Partial ROC is calculated following Peterson et al.
#' (2008; \url{http://dx.doi.org/10.1016/j.ecolmodel.2007.11.008}). This function is a modification
#' of the PartialROC funcion, available at \url{https://github.com/narayanibarve/ENMGadgets}.
#' @import Rcpp
#' @references Peterson, A.T. et al. (2008) Rethinking receiver operating characteristic analysis applications in ecological niche modeling. Ecol. Modell., 213, 63–72.
#' @examples
#' # Load a continuous model
#' conti_model <- raster::raster(system.file("extdata",
#'                                           "ambystoma_model.tif",
#'                                            package="ntbox"))
#' # Read validation (test) data
#' test_data <- read.csv(system.file("extdata",
#'                                   "ambystoma_validation.csv",
#'                                   package = "ntbox"))
#'
#' # Filter only presences as the Partial ROC only needs occurrence data
#' test_data <- dplyr::filter(test_data, presence_absence==1)
#' test_data <- test_data[,c("longitude","latitude")]
#'
#' partial_roc <- pROC(continuous_mod=conti_model,
#'                     test_data = test_data,
#'                     n_iter=1000,E_percent=5,
#'                     boost_percent=50,
#'                     parallel=FALSE)
#'
#' @importFrom purrr map_df
#' @import future
#' @useDynLib ntbox
#' @export


pROC <- function(continuous_mod,test_data,
                 n_iter=1000,E_percent=5,
                 boost_percent=50,
                 parallel=FALSE,ncores=4,rseed=FALSE){

  if (class(continuous_mod) == "RasterLayer") {
    if (continuous_mod@data@min == continuous_mod@data@max) {
      stop("\nModel with no variability.\n")
    }
    if (is.data.frame(test_data) || is.matrix(test_data)) {
      test_data <- stats::na.omit(raster::extract(continuous_mod,
                                                  test_data))

    }
    vals <- continuous_mod[!is.na(continuous_mod[])]
  }
  if(is.numeric(continuous_mod)){
    vals <- continuous_mod
    if (!is.numeric(test_data))
      stop("If continuous_mod is of class numeric,
           test_data must be numeric...")
  }
  ndigits <- proc_precision(mod_vals = vals,
                            test_data = test_data)

  tomult <- as.numeric(paste0("1e+",ndigits))
  test_value <- test_data*tomult
  test_value <- round(as.vector(test_value))

  vals2 <- round(vals*tomult)
  classpixels <- as.data.frame(base::table(vals2),
                               stringsAsFactors = F)
  names(classpixels) <- c("value",
                          "count")
  classpixels$value <- as.numeric(classpixels$value)
  classpixels <- data.frame(stats::na.omit(classpixels))
  value <- count <- totpixperclass <- NULL
  classpixels <- classpixels %>%
    dplyr::mutate(value  = rev(value),
                  count = rev(count),
                  totpixperclass = cumsum(count),
                  percentpixels = totpixperclass/sum(count)) %>%
    dplyr::arrange(value)

  #if(nrow(classpixels)>1500){
  #  classpixels <- classpixels %>%
  #    dplyr::sample_n(1500) %>% dplyr::arrange(value)
  #}

  error_sens <- 1 - (E_percent/100)
  models_thresholds <- classpixels[, "value"]
  fractional_area <- classpixels[, "percentpixels"]
  n_data <- length(test_value)
  n_samp <- ceiling((boost_percent/100) * (n_data))

  big_classpixels <- matrix(rep(models_thresholds,
                                each = n_samp),
                            ncol = length(models_thresholds))


  calc_aucDF <- function(big_classpixels,
                         fractional_area,
                         test_value,
                         n_data, n_samp,
                         error_sens,rseed=NULL) {
    if(is.numeric(rseed)) set.seed(rseed)

    rowsID <- sample(x = n_data,
                     size = n_samp,
                     replace = TRUE)

    test_value1 <- test_value[rowsID]
    omssion_matrix <- big_classpixels > test_value1
    sensibility <- 1 - colSums(omssion_matrix)/n_samp
    xyTable <- data.frame(fractional_area, sensibility)
    xyTable <- rbind(xyTable,c(0,0))
    xyTable <- xyTable[order(xyTable$fractional_area,
                             decreasing = F),]
    auc_model <- trapozoid_roc(xyTable$fractional_area,
                               xyTable$sensibility)
    if(error_sens>0){
      less_ID <- which(xyTable$sensibility <= error_sens)
      xyTable <- xyTable[-less_ID, ]
      auc_pmodel <- trapozoid_roc(xyTable$fractional_area,
                                  xyTable$sensibility)

      auc_prand <- trapozoid_roc(xyTable$fractional_area,
                                 xyTable$fractional_area)

    }
    else{
      auc_pmodel <- auc_model
      auc_prand <- 0.5

    }



    auc_ratio <- auc_pmodel/auc_prand
    auc_table <- data.frame(auc_model,
                            auc_pmodel,
                            auc_prand,
                            auc_ratio)
    return(auc_table)
  }

  if (parallel) {
    n_cores <- nc(ncores)

    furrr::furrr_options(packages = c("Rcpp","ntbox"))
    plan(multisession,workers=n_cores)
    options(future.globals.maxSize= 8500*1024^2)

    roc_env <- new.env()
    niter_big <- floor(n_iter/n_cores)
    n_runs <- rep(niter_big, n_cores)
    sum_n_runs <- sum(n_runs)
    n_runs[1] <- n_runs[1] + (n_iter - sum_n_runs)

    for (i in 1:length(n_runs)) {
      x <- as.character(i)
      roc_env[[x]] %<-% {
        library(Rcpp)
        x1 <- 1:n_runs[i]
        auc_matrix1 <- x1 %>%
          purrr::map_df(~calc_aucDF(big_classpixels,
                                    fractional_area,
                                    test_value,
                                    n_data, n_samp,
                                    error_sens,rseed=NULL))
      }
    }
    partial_AUC <- as.list(roc_env)
    rm(roc_env)
    partial_AUC <- do.call(rbind.data.frame,
                           partial_AUC)
    rownames(partial_AUC) <- NULL
    future::plan(future::sequential)
  }
  else {
    partial_AUC <- 1:n_iter %>%
      purrr::map_df(function(i){
        proc <- calc_aucDF(big_classpixels,
                           fractional_area,
                           test_value,
                           n_data,
                           n_samp,
                           error_sens,rseed = i)
      })
  }
  mauc <-  mean(partial_AUC$auc_model, na.rm = TRUE)
  maucp <- mean(partial_AUC$auc_ratio, na.rm = TRUE)
  proc <- sum(partial_AUC$auc_ratio <= 1, na.rm = TRUE)/
    length(partial_AUC$auc_ratio[!is.na(partial_AUC$auc_ratio)])

  p_roc <- c(mauc,maucp, proc)
  names(p_roc) <- c("Mean_AUC",
                    paste("Mean_pAUC_ratio_at_",
                          E_percent,
                          "%", sep = ""),
                    "P_value")
  p_roc_res <- list(pROC_summary = p_roc,
                    pROC_results = partial_AUC)
  return(p_roc_res)
}


proc_precision <- function(mod_vals,test_data){

  min_vals <- min(mod_vals,na.rm = TRUE)

  percentil_test <- unique(sort(stats::na.omit(test_data)))[2]


  #percentil_test <- stats::quantile(test_data,
  #                                  probs=0.1)
  partition_flag <- mean(c(min_vals,
                           percentil_test))
  fflag <- stringr::str_detect(partition_flag, "e")
  if (length(fflag)>0L && fflag) {
    ndigits <- stringr::str_split(partition_flag, "e-")[[1]]
    ndigits <- as.numeric(ndigits)[2] #- 1
  }
  else {
    med <- stringr::str_extract_all(partition_flag, pattern = "[0-9]|[.]")
    med <- unlist(med)
    med <- med[-(1:which(med == "."))]
    med1 <- which(med != 0)
    ndigits <- ifelse(med1[1] <= 2, 3, 4)
  }
  return(ndigits)
}

#' nc: Function to check the number of available cores
#'
#' @description nc is a helper function to check if the number of cores required to run a parallel process is less or equal to the total number of cores of the system.
#' @param ncores Number of cores for the parallel process.
#' @return  Returns an integer representing the number of cores that will be used to run a parallel process.
#' @details If ncores is bigger than the system's number of cores the function will return the system's number of cores. The functions that use this helper are \code{\link{mop}} and \code{\link{pROC}}.
#' @export
#' @examples
#' # Print the number of cores
#' print(nc(ncores=8))

nc<- function(ncores){
  nc_max <- future::availableCores()
  if(ncores > nc_max){
    warning(paste("The specified ncores are more than the system's ncores.",
                  "Runing with the system's ncores"))
    ncores <- nc_max
  }
  return(ncores)
}
