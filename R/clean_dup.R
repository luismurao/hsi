#' Function to clean duplicated longitude and latitude data
#' @description Clean duplicated longitude and latitude data using threshold distance
#'              which is a distance between points to be considered
#'              duplicates.
#' @param data A data.frame with longitude and latitude data
#' @param longitude A character vector of the column name of longitude.
#' @param latitude A character vector of the column name of latitude.
#' @param threshold A numerc value representig the eculidean distance between coordinates
#'            to be considered as a duplicate.
#' @return Returns a data.frame with coordinate data from species
#' @export
#' @examples
#' occ_rec <- read.csv(system.file("exdata/cardon_virtual.csv",
#'                                  package = "hsi"))
#' occ_clean <- hsi::clean_dup(occ_rec,longitude="x",
#'                             latitude="y",threshold=0.5)

clean_dup <- function(data,longitude,latitude,threshold=0.0){
  data <- data[!is.na(data[,longitude]),]
  dat_sp <- sp::SpatialPointsDataFrame(data[,c(longitude ,latitude)],data)
  dat_sp1 <- sp::remove.duplicates(dat_sp, zero = threshold)
  return(dat_sp1@data)
}


