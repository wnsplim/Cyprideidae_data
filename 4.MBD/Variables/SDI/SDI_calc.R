wd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)

library(chronosphere)
library(concaveman)
library(dplyr)
library(palaeoverse)
library(sf)
library(ggplot2)
library(patchwork)
library(scales)

occs <- read.table("./Cyprideidae_occ_standardized.txt", header = TRUE, sep = "\t")
Global_area <- read.table("./Global_area.txt", header = TRUE, sep = "\t")
NAfEu_area <- read.table("./NAfEu_area.txt", header = TRUE, sep = "\t")
WAfSAm_area <- read.table("./WAfSAm_area.txt", header = TRUE, sep = "\t")
  
maxage <- max(occs$Max_Ma)%/%5*5+5
minage <- min(occs$Min_Ma)%/%5*5
ages <- seq(minage, maxage, by = 5)
## Map for these ages are not available in PALEOMAP
ages <- setdiff(ages, c(100, 110, 115, 145))

CM_list <- list()
CS_list <- list()

for (i in 1:length(ages)) {
  cfiles <- list.files("./paleocoastlines_v7_shapefiles/CS/")
  sfiles <- list.files("./paleocoastlines_v7_shapefiles/CM/")
  CS_list[[i]] <- read_sf(dsn = "./paleocoastlines_v7_shapefiles/CS/", layer = gsub(".shp", "", cfiles[grep(paste0("^", ages[i], "Ma_CS_v7.shp$"), cfiles)]))
  CM_list[[i]] <- read_sf(dsn = "./paleocoastlines_v7_shapefiles/CM/", layer = gsub(".shp", "", sfiles[grep(paste0("^", ages[i], "Ma_CM_v7.shp$"), sfiles)]))
}

coastline_Global <- list()
land_Global<- list()
coastline_NAfEu <- list()
land_NAfEu<- list()
coastline_WAfSAm <- list()
land_WAfSAm <- list()

frame <- chronosphere::mapedge()

for(i in 1:length(ages)) {
  occs$age <- ages[i]
  Global_area$age <- ages[i]
  NAfEu_area$age <- ages[i]
  WAfSAm_area$age <- ages[i]
  
  rotation <- palaeorotate(occs, lng = "Longitude", lat = "Latitude", age = "age", model = "PALEOMAP") %>%
    dplyr::select(p_lat, p_lng)
  rotation_Global<- palaeorotate(Global_area, lng = "lng", lat = "lat", age = "age", model = "PALEOMAP") %>%
    dplyr::select(p_lat, p_lng)
  rotation_NAfEu<- palaeorotate(NAfEu_area, lng = "lng", lat = "lat", age = "age", model = "PALEOMAP") %>%
    dplyr::select(p_lat, p_lng)
  rotation_WAfSAm<- palaeorotate(WAfSAm_area, lng = "lng", lat = "lat", age = "age", model = "PALEOMAP") %>%
    dplyr::select(p_lat, p_lng)
  
  point_colors <- ifelse(occs$Original_Assignment %in% c("NAf", "Eu"), "navy",
                         ifelse(occs$Original_Assignment %in% c("WAf", "SAm"), "red",
                                ifelse(occs$Original_Assignment %in% c("SAm_Arg"), "orange",
                                       ifelse(occs$Original_Assignment %in% c("NAm"), "blue","grey"))))
  
  polygon_Global <- concaveman(st_as_sf(rotation_Global, coords = c("p_lng", "p_lat"), crs = 4326))
  polygon_NAfEu <- concaveman(st_as_sf(rotation_NAfEu, coords = c("p_lng", "p_lat"), crs = 4326))
  polygon_WAfSAm <- concaveman(st_as_sf(rotation_WAfSAm, coords = c("p_lng", "p_lat"), crs = 4326))
  
  intersection_Global <- st_intersection(polygon_Global, CS_list[[i]])
  intersection_NAfEu <- st_intersection(polygon_NAfEu, CS_list[[i]])
  intersection_WAfSAm <- st_intersection(polygon_WAfSAm, CS_list[[i]])
  
  # in km
  coastline_Global[[i]] <- sum(st_length(st_boundary(intersection_Global)))/1000
  coastline_NAfEu[[i]] <- sum(st_length(st_boundary(intersection_NAfEu)))/1000
  coastline_WAfSAm[[i]] <- sum(st_length(st_boundary(intersection_WAfSAm)))/1000

  # in km^2
  land_Global[[i]] <- sum(st_area(intersection_Global))/1e6
  land_NAfEu[[i]] <- sum(st_area(intersection_NAfEu))/1e6
  land_WAfSAm[[i]] <- sum(st_area(intersection_WAfSAm))/1e6
  
  plot(frame, col = "aliceblue")
  plot(CM_list[[i]], col = "lightsteelblue1", border = "lightsteelblue2", add = TRUE)
  plot(CS_list[[i]], col = "seashell", border = "steelblue2", add = TRUE)
  plot(frame, add = TRUE)
  
  plot(intersection_Global, col = "skyblue", border = "skyblue", add = TRUE)
  plot(intersection_NAfEu, col = "blue", border = "blue", add = TRUE)
  plot(intersection_WAfSAm, col = "pink", border = "pink", add = TRUE)
  
  points(rotation$p_lng, rotation$p_lat, pch = ".", cex = 3, col = point_colors)
}

coastline_Global <- as.numeric(gsub("\\[m\\]", "", unlist(coastline_Global)))
coastline_NAfEu <- as.numeric(gsub("\\[m\\]", "", unlist(coastline_NAfEu)))
coastline_WAfSAm <- as.numeric(gsub("\\[m\\]", "", unlist(coastline_WAfSAm)))

land_Global <- unlist(land_Global)
land_NAfEu <- unlist(land_NAfEu)
land_WAfSAm <- unlist(land_WAfSAm)

## Perform regression to estimate d value
df_Global <- data.frame(land = land_Global, coastline = coastline_Global) %>% 
  mutate(log_land = log10(land), log_coastline = log10(coastline))
df_NAfEu <- data.frame(land = land_NAfEu, coastline = coastline_NAfEu) %>% 
  mutate(log_land = log10(land), log_coastline = log10(coastline))
df_WAfSAm <- data.frame(land = land_WAfSAm, coastline = coastline_WAfSAm) %>% 
  mutate(log_land = log10(land), log_coastline = log10(coastline))

get_slope_ci <- function(df) {
  model <- lm(log_coastline ~ log_land, data = df)
  slope <- coef(model)[2]
  ci <- confint(model, level = 0.95)[2, ]
  return(sprintf("Slope: %.2f\n95%% CI: %.2f – %.2f", slope, ci[1], ci[2]))
}

plot_lm <- function(df, title, xpos, ypos, hjust) {
  slope_text <- get_slope_ci(df)
  
  ggplot(df, aes(x = log_land, y = log_coastline)) +
    geom_point(color = "blue", alpha = 0.6) +
    geom_smooth(method = "lm", color = "black", se = TRUE) +
    annotate("text", 
             x = xpos, 
             y = ypos,
             label = slope_text,
             hjust = hjust,
             vjust = 1,
             size = 10, fontface = "bold") +
    labs(title = title, x = "Land Area (log10)", y = "Coastline Length (log10)") +
    theme_bw(base_size = 24) +
    theme(panel.grid = element_blank(), plot.title = element_text(face = "bold", size = 30),
          axis.text = element_blank(), axis.ticks = element_blank())
}

plot_Global <- plot_lm(df_Global, "Global", max(df_Global$log_land), max(df_Global$log_coastline), 1)
plot_NAfEu <- plot_lm(df_NAfEu, "NAfEu", min(df_NAfEu$log_land), max(df_NAfEu$log_coastline), 0)
plot_WAfSAm <- plot_lm(df_WAfSAm, "WAfSAm", max(df_WAfSAm$log_land), max(df_WAfSAm$log_coastline), 1)

final_plot <- plot_Global + plot_NAfEu + plot_WAfSAm + plot_layout(ncol = 3)

print(final_plot)


SDI_Global <- coastline_Global/(2*sqrt(pi*land_Global))
SDI_NAfEu <- coastline_NAfEu/(2*sqrt(pi*land_NAfEu))
SDI_WAfSAm <- coastline_WAfSAm/(2*sqrt(pi*land_WAfSAm))


plot(ages, SDI_Global, type = "b", col = "blue", pch = 19, 
     xlab = "Age (Ma)", ylab = "Shoreline Development Index (SDI)", 
     xlim = rev(range(ages)), ylim = c(0, max(c(SDI_Global, SDI_NAfEu, SDI_WAfSAm))))

lines(ages, SDI_NAfEu, type = "b", col = "red", pch = 19)
lines(ages, SDI_WAfSAm, type = "b", col = "green", pch = 19)

legend("bottomright", legend = c("Global", "NAfEu", "WAfSAm"), 
       col = c("blue", "red", "green"), pch = 19, cex = 0.8)


SDI_Global_df <- data.frame(age = ages, SDI_Global = SDI_Global) %>% arrange(ages)
SDI_NAfEu_df <- data.frame(age = ages, SDI_NAfEu = SDI_NAfEu) %>% arrange(ages)
SDI_WAfSAm_df <- data.frame(age = ages, SDI_WAfSAm = SDI_WAfSAm) %>% arrange(ages)

write.table(SDI_Global_df, "SDI_Global.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(SDI_NAfEu_df, "SDI_NAfEu.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(SDI_WAfSAm_df, "SDI_WAfSAm.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

