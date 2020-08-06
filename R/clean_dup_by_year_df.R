#' clean_dup_by_year_df: Function to clean duplicated ocurrence data by year.
#' @description Clean duplicated occurence data by year using a threshold distance.
#' @param df_coords_year, A data.frame with coordinate records by year.
#' @param longitude A character vector of the column name of longitude.
#' @param latitude A character vector of the column name of latitude.
#' @param year_var  The name of the year variable.
#' @param threshold A numerc value representig the euclidean distance between coordinates.
#' @seealso  Use this version of the function \code{\link[hsi]{clean_dup_by_year}} for data.frames
#' @export
#' @examples
#' puma_recs <- read.csv(system.file("exdata/puma_coords_year.csv",
#'                                  package = "hsi"))
#' puma_clean_year <- clean_dup_by_year_df(df_coords_year = puma_recs,
#'                                      longitude = "longitude",
#'                                      latitude = "latitude",
#'                                      year_var = "year",
#'                                      threshold = 0.05)

clean_dup_by_year_df <- function(df_coords_year,longitude, latitude, year_var,threshold){

  df_coords_year[,year_var] <- as.numeric(as.character(df_coords_year[,year_var]))
  clean_by_year <- df_coords_year %>%
    split(df_coords_year[,year_var],drop=T) %>%
    purrr::map_df(~clean_dup(data = .x,
                             longitude = longitude,
                             latitude = latitude,
                             threshold = threshold))

  return(clean_by_year)
}
. <- NULL

