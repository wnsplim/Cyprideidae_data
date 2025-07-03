wd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)

library(dplyr)
library(palaeoverse)
library(readr)
library(terra)
library(sf)
library(exactextractr)
library(concaveman)

Global_area <- read.table("./Global_area.txt", header = TRUE, sep = "\t")
NAfEu_area <- read.table("./NAfEu_area.txt", header = TRUE, sep = "\t")
WAfSAm_area <- read.table("./WAfSAm_area.txt", header = TRUE, sep = "\t")

# One must download maps of Salles et al. (2023) corresponding to these ages from: http://www.hydroshare.org/resource/bbf65a3ba448456b8a7f2e2bb5e59733
ages <- c(81, 87, 92, 97, 103, 107, 111, 116, 122, 127, 131, 136, 142, 145, 149, 155, 160, 165, 168, 172)

init_df <- function() data.frame(age = numeric(0), value = numeric(0))
results <- list(
  LA_Global = init_df(), LA_NAfEu = init_df(), LA_WAfSAm = init_df(),
  PI_Global = init_df(), PI_NAfEu = init_df(), PI_WAfSAm = init_df(),
  BN_Global = init_df(), BN_NAfEu = init_df(), BN_WAfSAm = init_df(),
  BS_Global = init_df(), BS_NAfEu = init_df(), BS_WAfSAm = init_df()
)

for (age in ages) {
  
  Global_area$age <- age
  NAfEu_area$age <- age
  WAfSAm_area$age <- age
  
  get_polygon <- function(area_df) {
    palaeorotate(area_df, lng = "lng", lat = "lat", age = "age", model = "PALEOMAP") %>%
      dplyr::select(p_lat, p_lng) %>%
      st_as_sf(coords = c("p_lng", "p_lat"), crs = 4326) %>%
      concaveman()
  }
  
  polygons <- list(
    Global = get_polygon(Global_area),
    NAfEu = get_polygon(NAfEu_area),
    WAfSAm = get_polygon(WAfSAm_area)
  )

  lake_r   <- rast(nc_path, subds = "lakes")
  phydiv_r <- rast(nc_path, subds = "phydiv")
  basin_r  <- rast(nc_path, subds = "basin")
  
  lake_area_r <- lake_r %>%  
    cellSize(unit = "m") %>%
    mask(lake_r)
  
  basin_area_df <- lake_r %>%
    cellSize(unit = "m") %>%
    terra::zonal(basin_r, fun = "sum") %>%
    as.data.frame() %>%
    setNames(c("basin_id", "basin_area"))
  
  message(paste("Age:", age, "Ma"))
  
  for (region in names(polygons)) {  
    hull <- polygons[[region]]  
    
    message(paste("Region:", region))
    
    # Lakes
    lake_sum <- exact_extract(lake_area_r, hull, "sum")[[1]] / 1e+6
    results[[paste0("LA_", region)]] <- results[[paste0("LA_", region)]] %>%  
      bind_rows(data.frame(age = age, value = lake_sum))  
    message(paste("Lake area for", region, "-", round(lake_sum, 4), "km²"))
    
    # Physiography Index
    mean_phydiv <- exact_extract(phydiv_r, hull, "mean")[[1]]  
    results[[paste0("PI_", region)]] <- results[[paste0("PI_", region)]] %>%  
      bind_rows(data.frame(age = age, value = mean_phydiv))  
    message(paste("Physiography index for", region, "-", round(mean_phydiv, 4)))
    
    # Basins
    ex_ids <- terra::extract(basin_r, hull, cells = TRUE, ID = FALSE)[,1]  
    basin_ids <- unique(na.omit(ex_ids))  
    n_basins <- length(basin_ids)  
    mean_basin_area <- basin_area_df %>%  
        filter(basin_id %in% basin_ids) %>%  
        pull(basin_area) %>%  
        mean() / 1e+6
    
    results[[paste0("BN_", region)]] <- results[[paste0("BN_", region)]] %>%  
      bind_rows(data.frame(age = age, value = n_basins))  
    results[[paste0("BS_", region)]] <- results[[paste0("BS_", region)]] %>%  
      bind_rows(data.frame(age = age, value = mean_basin_area))  
    
    message(paste("Number of basins for", region, "-", n_basins))
    message(paste("Mean basin area for", region, "-", round(mean_basin_area, 4), "km²"))
  }  
  message(paste("Done:", age, "Ma"))
}


for (name in names(results)) {
  var_area <- strsplit(name, "_")[[1]]
  var <- var_area[1]
  area <- var_area[2]
  df <- results[[name]]
  colnames(df) <- c("age", paste0(var, "_", area))
  write.table(df, file = paste0(var, "_", area, ".txt"), sep = "\t", row.names = FALSE, quote = FALSE)
}

## Plots
# Lake Area (LA)
lake_area_Global <- results$LA_Global$value
lake_area_NAfEu <- results$LA_NAfEu$value
lake_area_WAfSAm <- results$LA_WAfSAm$value

plot(ages, lake_area_Global, type = "b", col = "blue", pch = 19,
     xlab = "Age (Ma)", ylab = "Lake Area (km²)", 
     xlim = rev(range(ages)), ylim = c(min(c(lake_area_Global, lake_area_NAfEu, lake_area_WAfSAm)) - 1000, max(c(lake_area_Global, lake_area_NAfEu, lake_area_WAfSAm)) + 1000))

lines(ages, lake_area_NAfEu, type = "b", col = "red", pch = 19)
lines(ages, lake_area_WAfSAm, type = "b", col = "green", pch = 19)

legend("bottomright", legend = c("Global", "NAfEu", "WAfSAm"), 
       col = c("blue", "red", "green"), pch = 19, cex = 0.8)


# Physiography Index (PI)
phydiv_Global <- results$PI_Global$value
phydiv_NAfEu <- results$PI_NAfEu$value
phydiv_WAfSAm <- results$PI_WAfSAm$value

plot(ages, phydiv_Global, type = "b", col = "blue", pch = 19,
     xlab = "Age (Ma)", ylab = "Physiography Index", 
     xlim = rev(range(ages)), ylim = c(0, 1))

lines(ages, phydiv_NAfEu, type = "b", col = "red", pch = 19)
lines(ages, phydiv_WAfSAm, type = "b", col = "green", pch = 19)

legend("bottomright", legend = c("Global", "NAfEu", "WAfSAm"), 
       col = c("blue", "red", "green"), pch = 19, cex = 0.8)


# Basin Count (BN)
basin_count_Global <- results$BN_Global$value
basin_count_NAfEu <- results$BN_NAfEu$value
basin_count_WAfSAm <- results$BN_WAfSAm$value

plot(ages, basin_count_Global, type = "b", col = "blue", pch = 19,
     xlab = "Age (Ma)", ylab = "Basin Count", 
     xlim = rev(range(ages)), ylim = c(min(c(basin_count_Global, basin_count_NAfEu, basin_count_WAfSAm)) - 1000, max(c(basin_count_Global, basin_count_NAfEu, basin_count_WAfSAm)) + 1000))

lines(ages, basin_count_NAfEu, type = "b", col = "red", pch = 19)
lines(ages, basin_count_WAfSAm, type = "b", col = "green", pch = 19)

legend("bottomright", legend = c("Global", "NAfEu", "WAfSAm"), 
       col = c("blue", "red", "green"), pch = 19, cex = 0.8)


# Basin Size (BS)
basin_size_Global <- results$BS_Global$value
basin_size_NAfEu <- results$BS_NAfEu$value
basin_size_WAfSAm <- results$BS_WAfSAm$value

plot(ages, basin_size_Global, type = "b", col = "blue", pch = 19,
     xlab = "Age (Ma)", ylab = "Mean Basin Size", 
     xlim = rev(range(ages)), ylim = c(min(c(basin_size_Global, basin_size_NAfEu, basin_size_WAfSAm), na.rm = TRUE) - 1000, max(c(basin_size_Global, basin_size_NAfEu, basin_size_WAfSAm), na.rm = TRUE) + 1000))

lines(ages, basin_size_NAfEu, type = "b", col = "red", pch = 19)
lines(ages, basin_size_WAfSAm, type = "b", col = "green", pch = 19)

legend("bottomright", legend = c("Global", "NAfEu", "WAfSAm"), 
       col = c("blue", "red", "green"), pch = 19, cex = 0.8)
