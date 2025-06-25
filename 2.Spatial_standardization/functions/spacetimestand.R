#' spacetimestand
#'
#' Function to spatially standardise occurrence data, with
#' the option to standardise through time. Two standardisation
#' metrics are implemented: longitude-latitude range and
#' minimum spanning tree length. Minimally, the function
#' requires longitude-latitude data to which the standardisation
#' routines will be applied. The function can also return spatial
#' summary statistics for the data - setting both cleaning routines
#' to FALSE therefore provides a useful way to calculate these
#' statistics for a given dataset without any standardisation.
#' Spatial metrics returned are the length of the minimum spanning
#' tree (mst_l) and the maximum great circle distance (gcm_mx) in
#' that tree in kilometres, the convex hull area (cha) and
#' convex hull perimeter length (chl) in square kilometres and
#' kilometres respectively, and the minimum (mn), maximum (mx),
#' mean (mn), median (md), upper (q1) and lower (q2) quartiles
#' and ranges of the latitudes (lat_r) and longitudes (lng_r)
#' in decimal degrees. If maximum and minimum ages are present,
#' the user can provide a vector of interval boundaries which
#' will be used to time-bin the data, and the standardisation
#' routines and/or summary calculations will be applied to each bin.
#' If chosen, the function can then calculate sampling-corrected
#' diversity to the spatially-standardised data using
#' coverage-based rarefaction, implemented by the estimateD
#' function of the iNEXT package.
#' @param x A dataframe containing, minimally, a column of
#' longitudes and a column of latitudes in decimal degrees
#' @param lng The name of the longitude column in x
#' @param lat The name of the latitude column in x
#' @param mst A logical determining if minimum spanning tree
#' length standardisation should be applied to the data
#' @param lnglat A logical determining if longitude-latitude
#' range standardisation should be applied to the data
#' @param stats A logical determining if spatial summary
#' statistics for the data should be returned
#' @param shapes A logical determining if the convex hulls
#' and minimum spanning trees for the data should be returned.
#' This is a convenience argument to provide Spatial* objects
#' useful for plotting
#' @param mst_l A numeric giving the target length of the
#' minimum spanning tree in thousands of kilometres after
#' standardisation
#' @param lng_r A numeric giving the target longitude range
#' in decimal degrees after standardisation. If NULL, the
#' function will internally assign a value equivalent to the
#' global longitudinal range i.e. no standardisation
#' @param lat_r A numeric giving the target latitude range
#' in decimal degrees after standardisation. If NULL, the
#' function will internally assign a value equivalent to the
#' global latitudinal range i.e. no standardisation
#' @param lng_tol The tolerance on the target longitude range
#' in decimal degrees
#' @param lat_tol The tolerance on the target latitude range
#' in decimal degrees. Set to equal to the lng tolerance by default
#' @param max_ma If not NULL, the column name in x of the
#' occurrence maximum ages, given in millions of years before
#' present
#' @param min_ma If not NULL, the column name in x of the
#' occurrence minimum ages, given in millions of years before
#' present
#' @param intervals If not NULL, a vector of time bin boundaries
#' given in millions of years before present, which will be used
#' to bin the occurrence data. If supplied, any returned data is
#' ordered from the oldest to the youngest bin
#' @param mode The desired binning method. One of 'intersect' or
#' 'exact'. If 'intersect' any occurrence ages which overlap with
#' a given interval will be taken. If 'exact' only occurrences
#' whose ages fall exactly in a given bin will be taken
#' @param area If not NULL, the area(s) within which data will
#' be subsampled before standardisation. This allows data from
#' specific regions to be standardised, rather than standardising
#' globally. This can be supplied as a two column dataframe/matrix
#' of longitude-latitude coordinates or a spatial polygon within
#' which the data will be subsampled, or a list of the same if
#' intervals if specified. Suitable windows can be generated
#' using @seealso spacetimewind
#' @param hex A vector of positive integers specifying a hexagonal
#' tesselation vector called by icosa::hexagrid. To reduce
#' computational burden, the supplied data is binned using this
#' hexagonal grid, then MST standardisation applied to the hexagon
#' centroids. As such, a igher tesselation will increase the
#' precision of the standardisation, but also the computation time
#' @param div A logical determining if sampling-corrected diversity
#' should be calculated. Note that this procedure can take a while
#' @param tax The name of the taxonomic names column in x
#' @param q A vector of coverage proportions, expressed as fractions
#' (0 < q < 1). Multiple proportions can be supplied
#' @param verbose A logical determining if function progress
#' should be reported
#' @return A list where the first element is the standardised data
#' (a dataframe subsampled from x, or a list of the same if intervals
#' was specified). Subsequent list elements are a dataframe of
#' spatial summary statistics if stats = TRUE, and a SpatialLines
#' object of minimum spanning trees and a SpatialPolygons object of
#' convex hulls if shapes = TRUE. If intervals were supplied, the
#' Spatial* objects will contain lines and polygons corresponding
#' to each bin
#' @export
#' @import sp sf geosphere data.table
#' @importFrom icosa hexagrid locate newsp centers
#' @importFrom igraph get.edgelist graph.adjacency as.directed as.undirected degree neighborhood
#' @importFrom GeoRange CHullAreaEarth
#' @importFrom stats na.omit as.dist
#' @importFrom ape mst
#' @importFrom grDevices chull
#' @importFrom iNEXT estimateD
#' @importFrom graphics hist


spacetimestand <- function(x, lng = "lng", lat = "lat", mst = TRUE, lnglat = TRUE, stats = TRUE, shapes = TRUE,
                           mst_l = NULL, lng_r = NULL, lat_r = NULL, lng_tol = 2, lat_tol = lng_tol,
                           max_ma = NULL, min_ma = NULL, intervals = NULL, mode = "intersect", area = NULL,
                           hex = c(4, 4), div = FALSE, tax = NULL, q = 0.5, verbose = TRUE) {

  library(sp)
  library(sf)
  library(data.table)
  library(geosphere)
  library(icosa)
  library(igraph)
  library(GeoRange)
  library(stats)
  library(ape)
  library(iNEXT)
  
# # 
# x = occs
# lng_r = lng[8]
# lat_r = lat[8]
# lng = "p_lng"
# lat = "p_lat"
# mst = TRUE
# lnglat = alnglat[8]
# stats = TRUE
# shapes = TRUE
# mst_l = amst[8]
# lng_tol = 2
# lat_tol = lng_tol
# max_ma = "max_ma"
# min_ma = "min_ma"
# intervals = stg
# mode = "intersect"
# area = alist[[8]]
# hex = c(4, 4)
# div = T
# tax = "genus"
# q = 0.5
# verbose = TRUE

  ######## ARG CHECKS ########

  # check that the requisite data minimally exists
  if(!exists("x")) {
    stop("Please supply x as a dataframe containing, minimally, a column of longitudes and a column of latitudes")
  }
  if(!is.data.frame(x) | ncol(x) < 2) {
    stop("Please supply x as a dataframe containing, minimally, a column of longitudes and a column of latitudes")
  }
  if(!is.null(intervals) & ncol(x) < 4) {
    stop("If intervals have been supplied, x must contain, minimally, a column of longitudes, a column of latitudes, a column of maximum ages and a column of minimum ages")
  }
  if(isTRUE(div) & ncol(x) < 3) {
    stop("If div is TRUE, x must contain, minimally, a column of longitudes, a column of latitudes, and a column of taxon names")
  }
  if(!is.null(intervals) & isTRUE(div) & ncol(x) < 5) {
    stop("If intervals have been supplied and div is TRUE, x must contain, minimally, a column of longitudes, a column of latitudes, a column of maximum ages, a column of minimum ages and a column of taxon names")
  }

  # check the coordinate data
  if(!is.character(lng) | length(lng) != 1) {
    stop("lng should be a character vector of length one, indicting the name of the longitude column in x")
  }
  if(!is.character(lat) | length(lat) != 1) {
    stop("lng should be a character vector of length one, indicting the name of the latitudes column in x")
  }
  if(!all(c(lng, lat) %in% colnames(x))) {
    stop("lng and lat should both be column names in x")
  }
  if(!all(is.numeric(x[,lng]), is.numeric(x[,lat]))) {
    stop("lng and lat must both refer to numeric columns in x")
  }
  if(any(c(is.na(x[,lng]), is.na(x[,lat])))) {
    stop("One or more longitude or latitude values are NA")
  }
  if(max(x[,lng]) > 180 | min(x[,lng]) < -180 | max(x[,lat] > 90 | min(x[,lat]) < -90)) {
    stop("Invalid coordinates present - all longitudes should be in range -180:180 and all latitudes in range -90:90")
  }

  # check methods
  if(!is.logical(mst) | length(mst) != 1) {
    stop("mst should be a logical of length 1, indicating whether MST length standardisation should be applied")
  }
  if(mst) {
    if(!is.numeric(mst_l) | length(mst_l) != 1) {
      stop("mst_l should be a numeric of length 1, giving the target minimum spanning tree length for the subsampled data in thousands of kilometres")
    }
    if(mst_l <= 0) {
      stop("mst_l must be positive")
    }
  }
  if(!is.logical(lnglat) | length(lnglat) != 1) {
    stop("lnglat should be a logical of length 1, indicating whether longitude-latitude range standardisation should be applied")
  }
  if(lnglat) {

    # set the global default if NULL
    if(is.null(lng_r)) {
      lng_r <- 359.9
    }
    if(!is.numeric(lng_r) | length(lng_r) != 1) {
      stop("lng_r should be a numeric of length 1, giving the target longitudinal range for the subsampled data in decimal degrees")
    }
    if(lng_r <= 0 | lng_r > 360) {
      stop("lng_r should be greater than 0 and less than 360")
    }
    if(!is.numeric(lng_tol) | length(lng_tol) != 1) {
      stop("lng_tol should be a numeric of length 1, giving the tolerance around the target longitudinal range for the subsampled data in decimal degrees")
    }
    if(lng_tol <= 0 | lng_tol > lng_r) {
      stop("lng_tol should be greater than 0 and less than lng_r")
    }
    # set the global default if NULL
    if(is.null(lat_r)) {
      lat_r <- 179.9
    }
    if(!is.numeric(lat_r) | length(lat_r) != 1) {
      stop("lat_r should be a numeric of length 1, giving the target latitudinal range for the subsampled data in decimal degrees")
    }
    if(lat_r <= 0 | lat_r > 180) {
      stop("lng_r should be greater than 0 and less than 180")
    }
    if(!is.numeric(lat_tol) | length(lat_tol) != 1) {
      stop("lat_tol should be a numeric of length 1, giving the tolerance around the target latitudinal range for the subsampled data in decimal degrees")
    }
    if(lat_tol <= 0 | lat_tol > lat_r) {
      stop("lat_tol should be greater than 0 and less than lat_r")
    }
  }

  # check time data if supplied
  if(!sum(c(is.null(max_ma), is.null(min_ma), is.null(intervals))) %in% c(0, 3)) {
    stop("If standardising the data through time, then max_ma, min_ma and intervals must all be supplied")
  }
  if(!is.null(max_ma)) {
    if(!is.character(max_ma) | length(max_ma) != 1) {
      stop("max_ma should be a character vector of length one, indicting the name of the maximum ages column in x")
    }
    if(!is.character(min_ma) | length(min_ma) != 1) {
      stop("min_ma should be a character vector of length one, indicting the name of the minimum ages column in x")
    }
    if(!all(c(max_ma, min_ma) %in% colnames(x))) {
      stop("lng and lat should both be column names in x")
    }
    if(!all(is.numeric(x[,max_ma]), is.numeric(x[,min_ma]))) {
      stop("max_ma and min_ma must both refer to numeric columns in x")
    }
    if(any(c(is.na(x[,max_ma]), is.na(x[,min_ma])))) {
      stop("One or more maximum or minimum ages is NA")
    }
    if(min(x[,max_ma]) < 0 | min(x[,min_ma]) < 0) {
      stop("One or more maximum ages is less than zero - ages should be expressed in million years before present")
    }
    if(!is.numeric(intervals)) {
      stop("Intervals must be supplied as a numeric vector of boundaries (minimally two unique elements), expressed in million years before present (positive values)")
    }
    if(any(is.na(intervals))) {
      stop("One or more elements of intervals is NA")
    }
    intervals <- unique(intervals)
    if(length(intervals) < 2 | min(intervals) < 0) {
      stop("Intervals must be supplied as a numeric vector of boundaries (minimally two unique elements), expressed in million years before present (positive values)")
    }
    #if(min(x[,min_ma]) < min(intervals) | max(x[,max_ma]) > max(intervals)) {
    #  if(verbose) {warning("Intervals does not cover the full time span of ages in x")}
    #}
    intervals <- intervals[order(intervals, decreasing = TRUE)]

    # as the default for intervals = NULL, create intervals as one single bin for the entire dataset
  } else {
    max_ma <- "use_max"
    min_ma <- "use_min"
    intervals <- c(2, 0)
  }

  # check area data if supplied
  if(!is.null(area)) {

    # check that the correct object type has been supplied
    if(!class(area) %in% c("data.frame", "matrix", "SpatialPolygons", "list")) {
      stop("Area should be a data.frame/matrix of longitude-latitude coordinates or a spatial polygon within which the data will be subsampled,
           or a list of the same if intervals if specified")
    }
    # coerce to list if needed
    if(is.data.frame(area) | is.matrix(area)) {
      area <- list(area)
    } else {
      area <- area
    }
    # ensure that the list is the same length as the number of intervals
    if(!length(area) %in% c(1, (length(intervals) - 1))) {
      stop("If intervals and area are both supplied, then area should be a data.frame/matrix of longitude-latitude coordinates or a spatial polygon
             within which the data will be subsampled in each time interval, or a list of the same with as many list elements as time intervals (i.e. length(intervals) - 1)")
    }
    if(length(area) == 1) {
      area_list <- lapply(1:(length(intervals) - 1), function(x) {x <- area})
      area <- area_list
    }

    # check validity of list items
    for(i in 1:length(area)) {

      ob <- area[[i]]
      if(is.data.frame(ob) | is.matrix(ob)) {
        ob <- as.matrix(unique(ob))
        if(ncol(ob) != 2) {
          stop("If supplying as a data.frame or matrix, area should consist of two columns denoting latitude and longitude respectively in decimal degrees")
        }
        if(nrow(ob) < 3) {
          stop("If supplying as a data.frame or matrix, area must contain at least three unique coordinates given in decimal degrees")
        }
        if(!is.numeric(ob[,1]) | !is.numeric(ob[,2])) {
          stop("Area must of class numeric")
        }
        if(any(c(is.na(ob[,1]), is.na(ob[,2])))) {
          stop("One or more values in area are NA")
        }
        if(max(ob[,1]) > 180 | min(ob[,1]) < -180 | max(ob[,2] > 90 | min(ob[,2]) < -90)) {
          stop("Invalid coordinates present in area - all longitudes should be in range -180:180 and all latitudes in range -90:90")
        }
        ob <- SpatialPolygons(list(Polygons(list(Polygon(ob)), ID = "A")), proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +lon_wrap=180"))
      }
      if(class(ob) == "SpatialPolygons") {
        if(length(ob@polygons) != 1) {
          stop("The SpatialPolygon should only contain a single polygon")
        }
        if(length(grep("\\+proj=longlat", as.character(ob@proj4string))) != 1) {
          stop("+proj=longlat not found in SpatialPolygon proj4string. Check that the SpatiaPolygon is in lng-lat projection")
        }
      }
      area[[i]] <- ob
    }

    # as the default for area = NULL, make a list the length of the bins, with each element as the surface of the entire earth
  } else {
    area <- lapply(1:(length(intervals) - 1), function(x) {SpatialPolygons(list(Polygons(list(Polygon(cbind(c(180, 180, -180, -180), c(90, -90, -90, 90)))), ID = as.character(x))), proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +lon_wrap=180"))})
  }

  # check diversity arguments if supplied
  if(div) {

    if(is.null(tax)) {
      stop("If div is TRUE, then tax must also be specified")
    }
    if(!is.character(tax) | length(tax) != 1) {
      stop("tax should be a character vector of length one, indicting the name of the taxon names column in x")
    }
    if(!tax %in% colnames(x)) {
      stop("tax should be a column name in x")
    }
    if(any(is.na(x[,tax]))) {
      stop("One or more values in column tax are NA")
    }
    if(any(is.na(q))) {
      stop("NA is not a permitted value in q")
    }
    if(!is.numeric(q) | any(q >= 1) | any(q <= 0)) {
      stop("q should be a vector of positive numerics (0 < q < 1) specifying the rarefaction values used in diversity estimation")
    }
  }

  # check minor arguments
  if(!is.numeric(hex)) {
    stop("hex should be a numeric vector of positive integers suitable for tesselation by icosa::hexagrid")
  }
  if(any(hex < 0 | !all(hex %% 1 == 0))) {
    stop("hex should be a numeric vector of positive integers suitable for tesselation by icosa::hexagrid")
  }
  if(!is.logical(stats) | length(stats) != 1) {
    stop("report should be a logical of length 1 determining if the spatial statistics of the subsamples should be returned or not")
  }
  if(!is.logical(shapes) | length(shapes) != 1) {
    stop("shapes should be a logical of length 1 determining if the minimum spanning trees and convex hulls of the subsamples should be returned or not")
  }
  if(!is.character(mode) | length(mode) != 1) {
    stop("mode should be a character of length one determining the time-binning method: one of 'intersect' or 'exact'")
  }
  if(!mode %in% c("intersect", "exact")) {
    stop("mode should be a character of length one determining the time-binning method: one of 'intersect' or 'exact'")
  }
  if(!is.logical(verbose) | length(verbose) != 1) {
    stop("verbose should be a logical of length 1 determining if the progress of the function should be reported in the console")
  }

  ######## INTERVAL LOOP ########

  # global variable workaround
  . <- NULL
  # copy to new object for the addition of new columns
  occs <- x
  # create hexagrid for computational reduction
  hx <- newsp(hexagrid(hex))
  fcentres <- icosa::centers(hx)
  pcentres <- icosa::locate(x = hx, y = occs[,c(lng, lat)], randomborder = TRUE)
  # add in row identifier to allow retrieval after division into bins
  occs$rnum <- 1:nrow(occs)
  # add in age defaults if unsupplied
  occs$use_max <- x$use_min <- 1
  # add in occurrence cells and cell centre coordinates
  occs$cell <- pcentres
  occs$clng <- fcentres[match(pcentres, rownames(fcentres)),1]
  occs$clat <- fcentres[match(pcentres, rownames(fcentres)),2]
  #occs$code <- paste(occs[,lng], occs[,lat], sep = "|")
  # convert to data.table to allow indexing to scale to very large datasets
  occs <- as.data.table(occs)

  # for each interval
  for(i in 1:(length(intervals) - 1)) {

    if(verbose) {
      if(i != 1) {cat("\r")}
      cat("Processing chunk ", i, "/", (length(intervals) - 1))
      if(i == (length(intervals) - 1)) {cat("\n")}
    }

    # initialise storage objects at start of loop
    if(i == 1) {
      win_stats <- win_pts <- win_mst <- win_ch <- win_div <- list()
      win_area <- area
    }


    # subset data temporally
    if(mode == "intersect") {
      o5 <- which(occs[[(min_ma)]] >= intervals[i])
      o6 <- which(occs[[(max_ma)]] <= intervals[i + 1])
      oa <- which(!(1:nrow(occs)) %in% unique(c(o5, o6)))
      pts <- as.data.frame(occs[oa,])
    }
    if(mode == "exact") {
      o5 <- which(pts[[(max_ma)]] <= intervals[i] & pts[[(min_ma)]] >= intervals[i + 1])
      pts <- as.data.frame(pts[oa,])
    }

    # subset data spatially
    poly <- win_area[[i]]
    pts <- pts[as.logical(point.in.polygon(pol.x = poly@polygons[[1]]@Polygons[[1]]@coords[,1], pol.y = poly@polygons[[1]]@Polygons[[1]]@coords[,2],
                          point.x = pts[,lng], point.y = pts[,lat])),]

    # get the cells covered by the data
    cpts <- pts[,c("clng", "clat")]
    cpts <- unique(cpts)

    ######## LNGLAT ########

    # do lnglat standardisation if called and if possible (i.e. above the 2 points minimumally needed to define a lnglat range)
    if(lnglat & nrow(pts) > 2) {

      # if the longitude range exceeds the threshold
      if(abs(diff(range(pts[,lng]))) > lng_r) {

        # tabulate the distribution of longitudes
        foo <- table(pts[,lng])

        # find range which is closest to the threshold, with the maximum data in that range
        sums <- list()
        for(j in 1:length(foo)) {
          upper <- which.min(abs((as.numeric(names(foo[j])) + lng_r) - as.numeric(names(foo[j:length(foo)]))))
          sms <- cumsum(foo[j:upper])
          sums[[j]] <- c(as.numeric(names(foo[c(j, upper)])), sms[length(sms)])
        }
        sums <- do.call(rbind, sums)

        # add in the difference between the target longitude range and the range of the solution
        sums <- cbind(sums, abs(lng_r - abs(sums[,1] - sums[,2])))

        # trim to solutions within the tolerance on acceptance
        sols <- which(sums[,4] < lng_tol)

        # if there are no solutions within the tolerance, take the closest solution
        if(length(sols) == 0) {sols <- which.min(sums[,4])}

        # get the occurrence totals within the solution ranges
        sums <- sums[sols, , drop = FALSE]

        # get the solution range with the most data
        lng_sol <- as.vector(sums[which.max(sums[,3]),1:2])
        lng_sol <- lng_sol[order(lng_sol)]

        # trim data to that range
        pts <- pts[pts[,lng] >= lng_sol[1] & pts[,lng] <= lng_sol[2],]
        cpts <- pts[,c("clng", "clat")]
        cpts <- unique(cpts)

      } else {
        lng_sol <- c(min(pts[,lng], max(pts[,lng])))
      }

      # if the latitude range exceeds the threshold
      if(abs(diff(range(pts[,lat]))) > lat_r) {

        # first subset to latitude points which fall outside the latitudes of the points defining the longitudes
        lng_def <- pts[which(pts[,lng] %in% lng_sol),]
        lat_u <- pts[which(pts[,lat] > max(lng_def[,lat])),]
        lat_l <- pts[which(pts[,lat] < min(lng_def[,lat])),]

        # if the limiting lats are also the limiting lngs, then bind the reduction points accordingly
        if(nrow(lat_u) == 0 & nrow(lat_l) == 0) {
          pts2 <- rbind(lng_def[which.min(lng_def[,lat]),],
                        lng_def[which.max(lng_def[,lat]),])
          lat_sol <- c(min(pts2[,lat]), max(pts2[,lat]))
        }
        if(nrow(lat_u) == 0 & nrow(lat_l) != 0) {
          pts2 <- rbind(lng_def[which.max(lng_def[,lat]),], lat_l)
          foo <- table(pts2[,lat])
          upper <- which.min(abs((as.numeric(names(foo[1])) + lat_r) - as.numeric(names(foo))))
          lat_sol <- as.numeric(names(foo[c(1, upper)]))
        }
        if(nrow(lat_u) != 0 & nrow(lat_l) == 0) {
          pts2 <- rbind(lng_def[which.min(lng_def[,lat]),], lat_u)
          foo <- table(pts2[,lat])
          upper <- which.min(abs((as.numeric(names(foo[1])) + lat_r) - as.numeric(names(foo))))
          lat_sol <- as.numeric(names(foo[c(1, upper)]))
        }
        if(nrow(lat_u) != 0 & nrow(lat_l) != 0) {
          pts2 <- rbind(lat_l, lat_u)
          foo <- table(pts2[,lat])
          sums <- list()
          for(j in 1:length(foo)) {
            upper <- which.min(abs((as.numeric(names(foo[j])) + lat_r) - as.numeric(names(foo[j:length(foo)]))))
            sms <- cumsum(foo[j:upper])
            sums[[j]] <- c(as.numeric(names(foo[c(j, upper)])), sms[length(sms)])
          }
          sums <- do.call(rbind, sums)
          sums <- cbind(sums, abs(lat_r - abs(sums[,1] - sums[,2])))
          # trim to solutions within a threshold acceptance
          sols <- which(sums[,4] < lat_tol[i])
          if(length(sols) == 0) {sols <- which.min(sums[,4])}
          sums <- sums[sols, , drop = FALSE]
          lat_sol <- as.vector(sums[which.max(sums[,3]),1:2])
          lat_sol <- lat_sol[order(lat_sol)]
        }

        # trim data to the range of the threshold solution
        pts <- pts[pts[,lat] >= lat_sol[1] & pts[,lat] <= lat_sol[2],]
        cpts <- pts[,c("clng", "clat")]
        cpts <- unique(cpts)
      }
    }

    ######## MST ########

    # do mst standardisation if called and if possible (2 cell points needed minimally to make a tree)
    if(mst & nrow(cpts) > 2) {

      # get initial mst length
      MSTd <- geosphere::distm(cpts[,1:2])
      MST <- get.edgelist(graph.adjacency(ape::mst(as.dist(MSTd))))
      len <- unlist(lapply(1:nrow(MST), function(x) {MSTd[as.numeric(MST[x,1]), as.numeric(MST[x,2])]}))
      len <- sum(len) / 1000 / 2

      # if the length is above the threshold
      if(len > mst_l) {

        # copy the cell points, get the hull points and remove
        sl <- cpts
        which_hull <- unique(c(which.max(sl[,1]), which.max(sl[,2]), which.min(sl[,1]), which.min(sl[,2])))
        hull <- sl[which_hull,]
        sl <- sl[-which_hull, , drop = FALSE]

        # tabulate occs per cell
        cls <- table(pts[,"cell"])

        # if there are non-hull points which can be removed
        if(nrow(sl) > 1) {

          # create loop breaker and point at beginning of loop
          keep_going <- TRUE
          iter <- 0

          # while loop to remove points iteratively
          while(len > mst_l & keep_going) {

            iter <- iter + 1
            # get the cells eligible for removal and their occurrence tallies
            loc_cell <- cls[match(locate(x = hx, y = sl), names(cls))]
            # get the adjacency matrix from the distance matrix, convert to edgelist (convert to undirected to allow redundant edge removal)
            gcd <- distm(sl) / 1000
            mst_loop <- as.directed(as.undirected(graph.adjacency(ape::mst(as.dist(gcd)))), mode = "arbitrary")
            # get the mst tips (graph degree = 1) with the node they link to, as an edgelist
            tips <- do.call(rbind.data.frame, neighborhood(mst_loop, nodes = as.vector(which(igraph::degree(mst_loop) == 1)), order = 1))
            # get the distances from the tips to their nearest points
            dists <- vector()
            for(k in 1:nrow(tips)) {
              dists[k] <- gcd[tips[k,1], tips[k,2]]
            }
            tips <- cbind.data.frame(tips, dists)
            colnames(tips) <- c("n1", "n2", "dist")
            # add in the occ counts for the tip cells and order
            tips$occs <- loc_cell[tips$n1]
            tips <- tips[order(tips$occs),]

            # subtract the length for the tip with the least data from the original mst length
            sol <- tips$n1[1]
            sol_d <- tips$dist[1]
            len_r <- len - sol_d

            # if the reduced length falls below the threshold, reassess length with the hull added back
            if(len_r <= mst_l) {

              # remove the tip and add back in the hull points
              sl2 <- sl[-sol, , drop = FALSE]
              sl2 <- rbind(sl2, hull)

              # recalculate the new mst length
              md <- geosphere::distm(sl2)
              ms <- get.edgelist(graph.adjacency(ape::mst(as.dist(md))))
              len2 <- unlist(lapply(1:nrow(ms), function(x) {md[as.numeric(ms[x,1]), as.numeric(ms[x,2])]}))
              len2 <- sum(len2) / 1000 / 2

              # if the reduced + hull solution falls below past the threshold
              if(len2 <= mst_l) {

                # and with a smaller distance to the threshold than the previous solution, take the new solution
                if(abs(mst_l - len2) < abs(mst_l - len)) {
                  sl <- sl[-sol, , drop = FALSE]
                  len <- len2

                  # otherwise take the previous solution and break the loop
                } else {
                  keep_going <- FALSE
                }

                # otherwise remove the points and proceed in the while loop
              } else {
                sl <- sl[-sol, , drop = FALSE]
                len <- len2
              }

              # otherwise remove the points and proceed in the while loop
            } else {
              sl <- sl[-sol, , drop = FALSE]
              len <- len_r
            }

            # if there are no more points available for removal, break the while loop
            if(nrow(sl) < 2) {
              keep_going <- FALSE
            }
          }

          # with the remaining points, add back the convex hull and recreate the mst
          sl_new <- rbind(sl, hull)
          # locate the cells of the old mst
          oldcells <- locate(x = hx, y = cpts[,c("clng", "clat")])
          # locate the cells of the new mst
          newcells <- locate(x = hx, y = sl_new)
          # subset pts to those in the new cells
          pts <- pts[pts$cell %in% newcells,]
          cpts <- pts[,c("clng", "clat")]
          cpts <- unique(cpts)
        }
      }
    }

    ######## STORE ########

    #win_pts[[i]] <- x[pts$rnum, ]
    #if(nrow(cpts) > 1) {
    #  MSTd <- geosphere::distm(cpts[, 1:2])
    #  MST <- get.edgelist(graph.adjacency(ape::mst(as.dist(MSTd))))
    #  begin <- cbind(cpts[as.numeric(MST[, 1]), 1], cpts[as.numeric(MST[,1]), 2])
    #  ends <- cbind(cpts[as.numeric(MST[, 2]), 1], cpts[as.numeric(MST[,2]), 2])
    #  mst_out <- vector("list", nrow(begin))
    #  mst_out <- SpatialLines(lapply(1:length(mst_out), function(x) {Lines(list(Line(rbind(begin[x, ], ends[x, ]))), as.character(x))}), proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +lon_wrap=180"))
    #  pts_ch <- as.matrix(cpts[chull(cpts[,1:2]), ])
    #  pts_ch <- rbind(pts_ch, pts_ch[1,])
    #  pts_ch <- SpatialPolygons(list(Polygons(list(Polygon(rbind(pts_ch, pts_ch[1,]))), ID = as.character(i))), proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +lon_wrap=180"))
    #} else {
    #  mst_out <- NULL
    #  pts_ch <- NULL
    #}
    #if(shapes & nrow(cpts) > 1) {
    #  wrap <- hist(cpts$clng, breaks = seq(from = -180, to = 180, by = 60), plot = FALSE)$counts
    #  if(all(wrap[2:5] == 0)) {
    #    mst_out <- recenter(mst_out) %>%
    #      as(., "sf") %>% st_wrap_dateline(options = c("WRAPDATELINE=YES")) %>%
    #      as_Spatial
    #    pts_ch <- recenter(pts_ch) %>%
    #      as(., "sf") %>% st_wrap_dateline(options = c("WRAPDATELINE=YES")) %>%
    #      as_Spatial
    #  }
    #  win_mst[[i]] <- mst_out
    #  win_ch[[i]] <- pts_ch
    #} else {
    #  win_mst[[i]] <- NULL
    #  win_ch[[i]] <- NULL
    #}

    win_pts[[i]] <- x[pts$rnum, ]
    if(nrow(cpts) > 1) {

      # hull
      pts_ch <- as.matrix(cpts[chull(cpts[,1:2]), ])
      pts_ch <- rbind(pts_ch, pts_ch[1,])
      pts_ch <- SpatialPolygons(list(Polygons(list(Polygon(rbind(pts_ch, pts_ch[1,]))), ID = as.character(i))), proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +lon_wrap=180"))

      # distance matrix for stats
      MSTd <- geosphere::distm(cpts[, 1:2])
      MST <- get.edgelist(graph.adjacency(ape::mst(as.dist(MSTd))))
      begin <- cbind(cpts[as.numeric(MST[, 1]), 1], cpts[as.numeric(MST[,1]), 2])
      ends <- cbind(cpts[as.numeric(MST[, 2]), 1], cpts[as.numeric(MST[,2]), 2])

      # create mst
      mst_out <- vector("list", nrow(begin))
      mst_out <- lapply(1:length(mst_out), function(x) {Lines(list(Line(rbind(begin[x, ], ends[x, ]))), as.character(x))})
      mst_out <- lapply(mst_out, function(x) {
        SpatialLines(list(x), proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +lon_wrap=180"))
      })
      # check each line for rewrapping
      to_wrap <- which(unlist(lapply(mst_out, function(x) {
        x@bbox[3] - x@bbox[1]
      })) > 180)
      # rewrap if needed
      if(length(to_wrap) > 0) {
        for(j in 1:length(to_wrap)) {
          ln <- recenter(mst_out[[to_wrap[j]]]) %>%
            as(., "sf") %>% st_wrap_dateline(options = c("WRAPDATELINE=YES")) %>% st_coordinates
          ln <- SpatialLines(list(Lines(lapply(list(ln[1:2,1:2], ln[3:4,1:2]), Line), ID = "A")),
                             proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +lon_wrap=180"))
          mst_out[[to_wrap[j]]] <- ln
        }
      }
      mst_out <- do.call(rbind, mst_out)
      # rewrapping test for the convex hull (approximated using the absence of low longitude points)
      wrap <- hist(cpts$clng, breaks = seq(from = -180, to = 180, by = 60), plot = FALSE)$counts
      if(all(wrap[2:5] == 0)) {
        pts_ch <- recenter(pts_ch) %>%
          as(., "sf") %>% st_wrap_dateline(options = c("WRAPDATELINE=YES")) %>%
          as_Spatial
      }
      # store if requested
      if(shapes) {
        win_mst[[i]] <- mst_out
        win_ch[[i]] <- pts_ch
      }

    } else {
      win_mst[[i]] <- NULL
      win_ch[[i]] <- NULL
    }

    ######## STATS ########

    if(stats) {

      # mst length (divide by two as lines are doubled up in edgelist)
      if(nrow(cpts) > 1) {
        mst_l_out <- unlist(lapply(1:nrow(MST), function(x) {MSTd[as.numeric(MST[x,1]), as.numeric(MST[x,2])]}))
        mst_l_out <- sum(mst_l_out) / 1000 / 2
        # max gcd length
        gcd <- max(na.omit(unlist(MSTd))) / 1000
      } else {
        mst_l_out <- gcd <- NA
      }
      if(nrow(cpts) > 2) {
        # convex hull area
        cha <- CHullAreaEarth(cpts[,1], cpts[,2])
        # convex hull perimeter length
        chl <- geosphere::perimeter(cpts[chull(cpts),]) / 1000
      } else {
        cha <- chl <- NA
      }
      if(nrow(win_pts[[i]]) > 1) {
        lt <- summary(win_pts[[i]][,lat])
        ltr <- diff(range(win_pts[[i]][,lat]))
        # lng stats
        ln <- summary(win_pts[[i]][,lng])
        lnr <- diff(range(win_pts[[i]][,lng]))
      } else {
        lt <- ln <- c(NA, NA, NA, NA, NA, NA)
        ltr <- lnr <- NA
      }
      # bind and return
      vals <- c(mst_l_out, gcd, cha, chl, lt, ln, ltr, lnr)
      win_stats[[i]] <- vals

      # finalise stats
      if(i == (length(intervals) - 1)) {
        win_stats <- do.call(rbind.data.frame, win_stats)
        names(win_stats) <- c("mst_l", "gcd_mx", "cha", "chl",
                              "lat_mn","lat_q1", "lat_md", "lat_mn", "lat_q2", "lat_mx", "lng_mn",
                              "lng_q1", "lng_md", "lng_mn", "lng_q2", "lng_mx", "lat_r", "lng_r")
      }
    }
  }

  ######## DIVERSITY ########

  if(div) {

    # get binned data
    region <- win_pts
    # raw diversity and data count
    n_div <- unlist(lapply(region, nrow))
    raw_div <- unlist(lapply(region, function(y) {length(unique(y[,tax]))}))

    # convert taxon counts to vector compatible with inext::estimateD
    for(j in 1:length(region)) {

      # if there is no data assign NA
      if(nrow(region[[j]]) == 0) {
        region[[j]] <- NA

        # otherwise tabulate the occurrences
      } else {
        gens <- table(region[[j]][,tax])
        gens <- gens[order(gens, decreasing = TRUE)]
        # if less than occs, assign NA as SQS will not work
        if(length(gens) < 3) {
          region[[j]] <- NA
          # otherwise
        } else {
          region[[j]] <- c(sum(gens), as.vector(gens))
        }
      }
    }
    # now record position of NA bins, then remove
    nas <- unlist(lapply(region, function(y) {
      is.na(y[1])
    }))
    if(sum(nas) != 0) {
      region <- region[!nas]
    }

    # now do sqs for the supplied quorum levels
    ests <- list()
    for(j in 1:length(q)) {
      
      if(length(region) == 0) {
        ests[[j]] <- data.frame(NA, NA, NA, NA)
        colnames(ests[[j]]) <- c(paste0("SC_", q[j]), paste0("qD_", q[j]), paste0("qD.LCL_", q[j]), paste0("warning_", q[j]))
        
      } else {
        estD <- suppressWarnings(estimateD(region, datatype = "incidence_freq", base = "coverage", level = q[j]))
        estD0 <- estD[estD$Order.q == 0,]
        # expand dataframe to include NA bins if needed
        if(sum(nas) != 0) {
          # take first row and blank it
          blank <- estD0[1, ,drop = FALSE]
          blank[] <- NA
          tmp <- blank[rep(1, length(nas)),]
          insert <- which(!nas)
          for(k in 1:length(insert)) {
            tmp[insert[k],] <- estD0[k,]
          }
          estD0 <- tmp
        }
        estD0$warning <- FALSE
        estD0$warning[which(estD0$t >= 3 * n_div)] <- TRUE
        # formatting
        estD0 <- estD0[5:8]
        colnames(estD0) <- paste0(colnames(estD0), paste0("_", q[j]))
        ests[[j]] <- estD0
        
      }
      if(verbose) {
        if(j != 1) {cat("\r")}
        cat(paste0(j, "/", length(q), " quorum levels done "))
        if(j == length(q)) {cat("\n")}
      }
    }
    
    win_div <- do.call(cbind, ests)
    # bind in sample size and raw diversity
    win_div$n <- n_div
    win_div$raw <- raw_div
  }

  ######## OUTPUT ########

  #outs <- c(TRUE, stats, shapes, shapes, shapes, div)
  #n_outs <- c("data", "stats", "mst", "chull", "window", "div")
  out <- list(win_pts, win_stats, win_mst, win_ch, win_area, win_div)
  names(out) <- c("data", "stats", "mst", "chull", "window", "div")
  out <- out[unlist(lapply(out, length)) > 0]
  return(out)
}
