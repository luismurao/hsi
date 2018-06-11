#' PartialROC - A modification of the function of the Barve & Barve (2016). ENMGadgets \url{https://github.com/narayanibarve}
#' Function PartialROC generates the area under the curve values using bootstrap method. PartialROC is a model evaluation tool, used for
#' continuous model outputs as compared to binary model outputs. This method is specially used for model trained using presence only data.
#' For more details refer DOI: 10.1016/j.ecolmodel.2007.11.008 and check ENMGadgets \url{https://github.com/narayanibarve}.
#' @param valData - Occurence validation data. Must have 3 columns SpName, Longitude, Latitude.
#' @param PredictionFile - It should be a raster class object of a continuous model output.
#' @param E - Amount of error admissible along the Y-axis, given the requirements and conditions of the study (by default =.05). Value should range between 0 - 1
#' @param RandomPercent - Occurrence points to be sampled randomly from the test data for bootstrapping.
#' @param NoOfIteration - Number of iteration for bootstrapping
#' @return OutputFile will have 4 columns, IterationNo, AUC_at_specified_value, AUC_AT_Random, AUC_Ratio. The first row will always have 0 th interation
#' which is the actual Area Under the Curve without bootstrapping. And the rest of the rows contains auc ratio for all the bootstrap.
#' @export


PartialROC <- function(valData, PredictionFile, E = 0.05,
                        RandomPercent, NoOfIteration)
{

  OmissionVal <- E
  #OutMat = matrix(0,nrow=NoOfIteration+1, ncol = 4)
  InRast = PredictionFile
  resacale_ras <- 10/raster::cellStats(InRast,max)
  ## Currently fixing the number of classes to 100. But later flexibility should be given in the parameter.
  InRast = round(InRast * resacale_ras)
  ClassPixels <- AreaPredictedPresence(InRast)
  Occur <- valData
  Occur = Occur[,-1]
  ExtRast = raster::extract(InRast, Occur)

  OccurTbl = cbind(Occur, ExtRast)
  OccurTbl = OccurTbl[which(is.na(OccurTbl[,3]) == FALSE),]

  PointID = seq(1:nrow(OccurTbl))
  OccurTbl = cbind(PointID, OccurTbl)
  names(OccurTbl)= c("PointID", "Longitude",
                     "Latitude", "ClassID")

  output_auc <- parallel::mclapply( 1:(NoOfIteration),
                                    function(x) auc_comp(x,OccurTbl,
                                                         RandomPercent,
                                                         OmissionVal,
                                                         ClassPixels))
  pRoc <- data.frame(t(sapply(output_auc,c)))
  return( pRoc)
}



#' Helper function to compute partial AUC, AUC ratio.
#' @param IterationNo number of iteration to compute partial AUC values
#' @param OccurTbl Validation data. Must have 3 columns SpName, Longitude, Latitude.
#' @param RandomPercent Occurrence points to be sampled randomly from the test data for bootstrapping.
#' @param OmissionVal 1-E.
#' @param ClassPixels Pixel classes.

auc_comp <- function(IterationNo,OccurTbl,RandomPercent,OmissionVal,ClassPixels){
  ClassID <- NULL
  n <- NULL
  OccuSumBelow <- NULL

  if (IterationNo > 0){
    ll = sample(nrow(OccurTbl),
                round(RandomPercent/100 * nrow(OccurTbl)),
                replace=TRUE)
    OccurTbl1 = OccurTbl[ll,]
  }
  else
    OccurTbl1 = OccurTbl
  OccurINClass <- OccurTbl1 %>% dplyr::group_by(ClassID) %>%
    dplyr::count() %>% dplyr::arrange(dplyr::desc(ClassID))
  OccurINClass <- OccurINClass %>%
    dplyr::ungroup() %>% dplyr::mutate(OccuSumBelow= cumsum(n)) %>%
    dplyr::mutate(Percent= OccuSumBelow/nrow(OccurTbl1))
  OccurINClass <- as.data.frame(OccurINClass)

  names(OccurINClass) = c("ClassID","OccuCount",
                          "OccuSumBelow", "Percent")

  XYTable = GenerateXYTableb(ClassPixels,OccurINClass)
  AreaRow = CalculateAUC(XYTable, OmissionVal, IterationNo)

  names(AreaRow) <- c( "IterNum","AUC_at_Value_0.95",
                       "AUC_at_0.5", "AUC_ratio")
  return(AreaRow)
}


#' Helper function to compute the area (number of pixels) that a certain threshold has.
#' @param InRast A raster class object of a continuous model output.

AreaPredictedPresence <- function(InRast)
{
  ### Now calculate proportionate area predicted under each suitability
  ClassPixels = raster::freq(InRast)
  ### Remove the NA pixels from the table.
  if (is.na(ClassPixels[dim(ClassPixels)[1],1])== TRUE)
  {
    ClassPixels = ClassPixels[-dim(ClassPixels)[1],]
  }

  ClassPixels = ClassPixels[order(nrow(ClassPixels):1),]
  TotPixPerClass = cumsum(ClassPixels[,2])
  PercentPixels = TotPixPerClass / sum(ClassPixels[,2])

  ClassPixels = cbind(ClassPixels, TotPixPerClass, PercentPixels)
  ClassPixels = ClassPixels[order(nrow(ClassPixels):1),]
  return(ClassPixels)
}


#' Helper function to compute the area (number of pixels) that a certain threshold has.
#' @param ClassPixels Pixel threshold class .
#' @param OccurINClass Ocurrence points that lies in a certain class.


GenerateXYTableb <- function(ClassPixels,OccurINClass){
  XYTable = data.frame(ClassPixels[,c(1,4)])
  names(XYTable) <- c("ClassID","PercentPixels")
  XYTable <-  suppressMessages(dplyr::full_join(XYTable,OccurINClass))

  XYTable$Percent[1] <- 1
  XYTable$Percent[is.na(XYTable$Percent)] <- OccurINClass[1,"Percent"]
  XYTable <- XYTable[,c("ClassID","PercentPixels","Percent")]
  XYTable <- rbind(XYTable,c(nrow(XYTable)+1,0,0))
  names(XYTable) = c("ClassID", "XCoor", "YCoor")
  return(XYTable)

}

#' Helper function to compute AUC (partialAUC, AUC at Random, AUC ratio) values
#' @param XYTable A table with the output of the function GenerateXYTableb
#' @param OmissionVal Omission value.
#' @param IterationNo Number of boostrap interation.

CalculateAUC <- function(XYTable, OmissionVal, IterationNo)
{
  ## if OmissionVal is 0, then calculate the complete area under the curve. Otherwise calculate only partial area


  if (OmissionVal > 0)
  {
    PartialXYTable = XYTable[which(XYTable[,3] >= OmissionVal),]
    ### Here calculate the X, Y coordinate for the parallel line to x-axis depending upon the OmissionVal
    ### Get the classid which is bigger than the last row of the XYTable and get the XCor and Ycor for that class
    ### So that slope of the line is calculated and then intersection point between line parallel to x-axis and passing through
    ### ommissionval on Y-axis is calculated.
    PrevXCor = XYTable[which(XYTable[,1]==PartialXYTable[nrow(PartialXYTable),1])+1,2]
    PrevYCor = XYTable[which(XYTable[,1]==PartialXYTable[nrow(PartialXYTable),1])+1,3]
    XCor1 = PartialXYTable[nrow(PartialXYTable),2]
    YCor1 = PartialXYTable[nrow(PartialXYTable),3]
    ## Calculate the point of intersection of line parallel to x-asiz and this line. Use the equation of line
    ## in point-slope form y1 = m(x1-x2)+y2
    Slope = (YCor1 - PrevYCor) / (XCor1 - PrevXCor)
    YCor0 = OmissionVal
    XCor0 = (YCor0 - PrevYCor + (Slope * PrevXCor)) / Slope
    ### Add this coordinate in the PartialXYTable with classid greater than highest class id in the table.
    ### Actually class-id is not that important now, only the place where we add this xcor0 and ycor0 is important.
    ### add this as last row in the table
    PartialXYTable = rbind(PartialXYTable,
                           c(PartialXYTable[nrow(PartialXYTable),1]+1,
                             XCor0, YCor0))
  }
  else
  {
    PartialXYTable = XYTable
  } ### if OmissionVal > 0

  ## Now calculate the area under the curve on this table.
  XCor1 = PartialXYTable[nrow(PartialXYTable),2]
  YCor1 = PartialXYTable[nrow(PartialXYTable),3]
  AUCValue = 0
  AUCValueAtRandom = 0
  for (i in (nrow(PartialXYTable)-1):1)
  {
    XCor2 = PartialXYTable[i,2]
    YCor2 = PartialXYTable[i,3]

    # This is calculating the AUCArea for 2 point trapezoid.
    TrapArea = (YCor1 * (abs(XCor2 - XCor1))) +
      (abs(YCor2 - YCor1) * abs(XCor2 - XCor1)) / 2
    AUCValue = AUCValue + TrapArea
    # now caluclate the area below 0.5 line.
    # Find the slope of line which goes to the point
    # Equation of line parallel to Y-axis is X=k and equation of line at 0.5 is y = x
    TrapAreaAtRandom = (XCor1 * (abs(XCor2 - XCor1))) +
      (abs(XCor2 - XCor1) * abs(XCor2 - XCor1)) / 2
    AUCValueAtRandom = AUCValueAtRandom + TrapAreaAtRandom
    XCor1 = XCor2
    YCor1 = YCor2

  }

  NewRow = c(IterationNo, AUCValue,
             AUCValueAtRandom,
             AUCValue/AUCValueAtRandom)

  return(NewRow)

}


