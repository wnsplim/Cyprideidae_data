#' spacetimewind
#'
#' Function to create a series of spatially constant
#' windows across a given number of time bins, with
#' the ability to iteratively offset the position of
#' the window from bin to bin in a (semi) regular
#' fashion. Can be used to create a series of windows
#' compatible with @seealso spacetimestand. In each
#' bin, the shift is multiplied by the bin index - 1,
#' resulting in a regular shift through time with no
#' movement in the first bin. The shift is defined
#' in the arguments in decimal degrees as these are
#' more intuitive to work with on a global scale. To
#' avoid distortion of distances at high longitudes
#' and latitudes, the windows are converted internally
#' from longitude-latitude coordinates to Robinson
#' projected coordinates, and the shift in decimal
#' degrees converted to metres using the approximation
#' of 1 degree = 111 km. The shift in metres is
#' applied to the projected coordinates, then the
#' coordinates converted back to longitude and latitude
#' before returning as a list of SpatialPolygons.
#' Regardless of the shift or multiplier, windows will
#' not shift beyond the boundaries of the poles.
#' Dateline wrapping is supported, however, for windows
#' that shift across the meridians
#' @param x A two column, numeric dataframe or matrix
#' of longitudes and latitudes used to define the
#' spatial window. This window will be returned as
#' the first in the series (i.e. with no shift)
#' @param shift A vector of length two defining the
#' longitude and latitude shift. Negative values for
#' longitude and latitude will produce westward and
#' southward movements respectively
#' @param bins A positive integer defining the number
#' of time bins over which the window will be shifted
#' @param lng_m If not NULL, a positive numeric which
#' will act as an additional multiplier during the bin
#' shift i.e. the longitude shift will be of the form
#' lng_shift \* bin_index \* lng_m, rather than just
#' lng_shift \* bin_index. This allows for a non-linear
#' shift in the window, although this parameter is
#' quite sensitive and more trial and error may be
#' needed to produce the desired series of windows
#' @param lat_m If not NULL, a positive numeric which
#' will act as an additional multiplier during the bin
#' shift i.e. the latitude shift will be of the form
#' lat_shift \* bin_index \* lat_m, rather than just
#' lat_shift \* bin_index. This allows for a non-linear
#' shift in the window, although this parameter is
#' quite sensitive and more trial and error may be
#' needed to produce the desired series of windows
#' @return A list containing as the spatial windows
#' for each bin. Each element is a SpatialPolygon in
#' unprojected coordinates
#' @export
#' @import sp sf
#' @importFrom methods as
#' @importFrom magrittr %>%
#'
spacetimewind <- function(x, shift, bins, lng_m = NULL, lat_m = NULL) {
  
  library(sp)
  library(sf)
  library(magrittr)

  # check coordinates
  if(!exists("x")) {
    stop("Please supply x as a two-column dataframe or matrix of longitudes and latitudes")
  }
  if(!class(x)[1] %in% c("data.frame", "matrix")) {
    stop("Please supply x as a two-column dataframe or matrix of longitudes and latitudes")
  }
  if(ncol(x) != 2) {
    stop("Please supply x as a two-column dataframe or matrix of longitudes and latitudes")
  }
  if(!is.numeric(x[,1]) | !is.numeric(x[,2])) {
    stop("x must contain numeric data only")
  }
  if(any(is.na(x))) {
    stop("One or more elements of x is NA")
  }
  if(max(x[,1]) > 180 | min(x[,1]) < -180 | max(x[,2] > 90 | min(x[,2]) < -90)) {
    stop("Invalid coordinates present - all longitudes should be in range -180:180 and all latitudes in range -90:90")
  }

  # check shift parameters
  if(!exists("shift")) {
    stop("Please supply shifts as a numeric vector of length two, containing the longitude and latitude shift to be applied to the coordinates in x")
  }
  if(!is.vector(shift)) {
    stop("Shifts is not a vector")
  }
  if(!is.numeric(shift)) {
    stop("Shifts is not numeric")
  }
  if(any(is.na(shift))) {
    stop("NA values are not allowed in shifts")
  }
  if(!exists("bins")) {
    stop("bins should be a single integer specifying the number of times the window should be moved by shift")
  }
  if(!is.numeric(bins) | bins < 1 | !bins %% 1 == 0) {
    stop("bins should be a single integer specifying the number of times the window should be moved by shift")
  }
  if(!is.null(lng_m)) {
    if(!is.numeric(lng_m) | length(lng_m) != 1 | lng_m < 0) {
      stop("lng_m should be a single positive numeric specifying the multiplier to be applied to longitude component of shift")
    }
    lng_m <- rep(lng_m, bins)
  } else {
    lng_m <- 1 / 1:bins
  }
  if(!is.null(lat_m)) {
    if(!is.numeric(lat_m) | length(lat_m) != 1 | lat_m < 0) {
      stop("lat_m should be a single positive numeric specifying the multiplier to be applied to latitude component of shift")
    }
    lat_m <- rep(lat_m, bins)
  } else {
    lat_m <- 1 / 1:bins
  }

  # global variable workaround
  . <- NULL

  # define equirectangular degrees approximately in metres (1d = 111 km = 111*10^3 m)
  shift <- shift * 111*(10^3)

  # define projections
  equi <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +no_defs")
  robin <- CRS("+proj=robin +datum=WGS84 +ellps=WGS84 +lon_0=0 +x_0=0 +y_0=0 +units=m +no_defs")

  # make spatial polygon, convert to sf, project to robinson
  poly <- SpatialPolygons(list(Polygons(list(Polygon(x)), ID = "A")), proj4string = equi) %>%
    st_as_sfc %>%
    st_transform(crs = robin)

  # do shift through time in metres of the form shift * iteration of the shift y * multiplier *_m
  poly_set <- lapply(1:bins, function(y) {poly + c((y - 1) * (shift[1] * (y * lng_m[y])), (y - 1) * (shift[2] * (y * lat_m[y])))})

  # for each shifted polygon, restore Robinson CRS, reproject to lng-lat, return to spatial polygon with dateline wrapping
  for (i in 1:length(poly_set)) {
    poly_set[[i]] <- `st_crs<-`(poly_set[[i]], robin) %>%
      st_transform(., crs = equi) %>% st_cast %>%
      lapply(., function(y) {Polygons(lapply(seq_along(y), function(z) {Polygon(y[[z]], z > 1)}), "ID")}) %>%
      SpatialPolygons(., proj4string = equi)

    # if self intersection is detected, assume this arises from dateline wrapping issues and correct
    if(suppressWarnings(!sf::st_is_valid(st_as_sf(poly_set[[i]])))) {
      poly_set[[i]] <- recenter(poly_set[[i]]) %>%
        as(., "sf") %>% st_wrap_dateline(options = c("WRAPDATELINE=YES")) %>%
        as_Spatial
    }
  }

  # return list of spatial polygons
  return(poly_set)
}
