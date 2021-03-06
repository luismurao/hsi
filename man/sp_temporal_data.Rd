% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sp_temporal_data.R
\name{sp_temporal_data}
\alias{sp_temporal_data}
\title{Function create a Species Temporal Data object (STD object).}
\usage{
sp_temporal_data(
  occs = NA,
  longitude = NULL,
  latitude = NULL,
  sp_year_var = NA,
  layers_by_year_dir = NA,
  layers_ext = "*.tif$",
  reclass_year_data = TRUE
)
}
\arguments{
\item{occs}{A data.frame or a SpatialPointsDataFrame with coordinates of the occurrence records. The data must have a time variable indicating the year of each observation.}

\item{longitude}{If occs is a data.frame the user must indicate the variable name of logitude data.}

\item{latitude}{If occs is a data.frame the user must indicate the variable name of latitude data.}

\item{sp_year_var}{A time variable indicating the year of each observation.}

\item{layers_by_year_dir}{El directorio que contiene las carpetas donde estan las capas raster que se usaran para modelar; estas deben estar clasificadas por anio.}

\item{layers_ext}{Extension de los archivos de las capas raster.}

\item{reclass_year_data}{Logical. If TRUE, occurrence data whose years are not in the environmental layers will be reclassified to the nearest year in the layers. Other wise they will not be used for modeling.}
}
\value{
Returns un  sp.temporal.modeling object (list) con las coordenadas de los puntos de presencia, los anios las observaciones, el path de las capas temporales.
}
\description{
This function creates Species Temporal Data object.
}
