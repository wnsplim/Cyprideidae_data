## This is a customization of a script "2_spatial_standardisation.R" from Flannery-Sutherland et al. (2024; Nat. Comm. 15:5382).
## Please consult the original script for details.

wd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)

library(chronosphere) # v 0.4.1
library(sp)
library(sf)
library(geosphere)
library(data.table)
library(icosa)
library(igraph)
library(GeoRange)
library(ape)
library(iNEXT)
library(zen4R)
library(readr)
library(raster)
source("functions/spacetimewind.R")
source("functions/spacetimestand.R")

occs <- read.table("data/Cyprideidae_occ.txt", header = TRUE, sep = "\t")

stg <- c(174.7, 161.5, 145.0, 139.8, 132.6, 125.77, 121.4, 113.0, 100.5, 89.8)
mpt <- stg[1:(length(stg) - 1)] - abs((diff(stg)) / 2)

fives <- 5 * round(mpt / 5)
## Map for these ages are not available in PALEOMAP
fives[which(fives == 115)] <- 120

shelf_list <- list()
coast_list <- list()

for(i in 1:length(fives)) {
  cfiles <- list.files("data/coast/")
  sfiles <- list.files("data/shelf/")
  coast_list[[i]] <- read_sf(dsn = "data/coast", layer = gsub(".shp", "", cfiles[grep(paste0("^", fives[i], "Ma_CS_v7.shp$"), cfiles)]))
  shelf_list[[i]] <- read_sf(dsn = "data/shelf", layer = gsub(".shp", "", sfiles[grep(paste0("^", fives[i], "Ma_CM_v7.shp$"), sfiles)]))
}

frame <- chronosphere::mapedge()

wind_NAfEu <- spacetimewind(cbind(c(-20,  -20,  45, 63, 40, 20),
                                c(19, 50, 50, 29, -5, -5)),
                          shift = c(0, 1.2), bins = length(stg) - 1)

wind_WAfSAm <- spacetimewind(cbind(c(-20, -20, 0, 30, 30, 0),
                                   c(-4, -40, -40, -30, -10, -4)),
                             shift = c(0, 1.2), bins = length(stg) - 1)

for(i in 1:length(shelf_list)) {
  binned <- occs[-(which(occs$Max_Ma <= stg[i + 1] | occs$Min_Ma >= stg[i])),]
  binned <- occs[which(occs$Max_Ma <= stg[i] & occs$Max_Ma >= stg[i + 1] | occs$Min_Ma <= stg[i] & occs$Min_Ma >= stg[i + 1]),]
  
  plot(frame, col = "aliceblue")
  plot(shelf_list[[i]], col = "lightsteelblue1", border = "lightsteelblue2", add = TRUE)
  plot(coast_list[[i]], col = "seashell", border = "steelblue2", add = TRUE)
  plot(frame, add = TRUE)
  points(binned$Paleo_Longitude, binned$Paleo_Latitude, pch = ".", cex = 3)
  
  plot(wind_NAfEu[[i]], border = "mediumseagreen", lwd = 2, add = TRUE)
  plot(wind_WAfSAm[[i]], border = "orange", lwd = 2, add = TRUE)
}

alist <- list(NULL, wind_NAfEu, wind_WAfSAm)
nms <-     c("Global", "NAfEu", "WAfSAm") 
amst <-    c(10000, 6000, 3000)
alnglat <- c(T, T, T)
lng <-     c(180, 90, 90)
lat <-     c(120, 90, 90)

amst_unstd <-    c(1000000000, 1000000000, 1000000000)
lng_unstd <-     c(360, 360, 360)
lat_unstd <-     c(180, 180, 180)

for(s in 1:length(alist)) {
  unstandardised <- spacetimestand(x = occs, lng = "Paleo_Longitude", lat = "Paleo_Latitude", area = alist[[s]],
                                   mst = FALSE, lnglat = FALSE, intervals = stg, max_ma = "Max_Ma", min_ma = "Min_Ma")
  
  standardised <- spacetimestand(x = occs, lng = "Paleo_Longitude", lat = "Paleo_Latitude", area = alist[[s]],
                                 mst = TRUE, mst_l = amst[s], lnglat = alnglat[s], lng_r = lng[s], lat_r = lat[s],
                                 intervals = stg, max_ma = "Max_Ma", min_ma = "Min_Ma", div = T, tax = "Taxon_name")
  
  
  cairo_pdf(paste0("plots/", nms[s], "_standardization.pdf"), width = 7, height = 5)
  
  layout(matrix(1:4, ncol = 2, byrow = T), heights = c(1, 1.25), widths = c(1, 1.2))
  par(mar = c(0.5, 4.5, 1, 0.5))
  
  if(all(is.na(unstandardised$stats$mst_l))) {
    plot(NA, NA, xlim = c(stg[1], stg[length(stg)]), xlab = "", ylab = "MST length (km)", ylim = c(0, 0), xaxt = "n")
  } else {
    plot(mpt, unstandardised$stats$mst_l, xlim = c(stg[1], stg[length(stg)]), xlab = "", ylab = "MST length (km)", xaxt = "n",
         ylim = c(min(na.omit(standardised$stats$mst_l)), max(na.omit(unstandardised$stats$mst_l))), type = "l")
    lines(mpt, standardised$stats$mst_l, xlim = c(stg[1], stg[length(stg)]), lty = 2)
  }
  axis(1, labels = F)
  
  par(mar = c(0.5, 4.5, 1, 4.5))
  
  plot(mpt, unstandardised$stats$lng_r, xlim = c(stg[1], stg[length(stg)]), xlab = "", ylab = "Longitude range (deg)", xaxt = "n",
       ylim = c(min(na.omit(standardised$stats$lng_r)), max(na.omit(unstandardised$stats$lng_r))), type = "l")
  par(new = TRUE)
  plot(mpt, standardised$stats$lng_r, xlim = c(stg[1], stg[length(stg)]), axes = FALSE, xlab = "", ylab = "", xaxt = "n", yaxt = "n",
       ylim = c(min(na.omit(standardised$stats$lng_r)), max(na.omit(unstandardised$stats$lng_r))), type = "l", lty = 2)
  par(new = TRUE)
  plot(mpt, unstandardised$stats$lat_r, xlim = c(stg[1], stg[length(stg)]), axes = FALSE, xlab = "", ylab = "", xaxt = "n", yaxt = "n",
       ylim = c(min(na.omit(standardised$stats$lat_r)), max(na.omit(unstandardised$stats$lat_r))), type = "l", col = 2)
  par(new = TRUE)
  plot(mpt, standardised$stats$lat_r, xlim = c(stg[1], stg[length(stg)]), axes = FALSE, xlab = "", ylab = "", xaxt = "n", yaxt = "n",
       ylim = c(min(na.omit(standardised$stats$lat_r)), max(na.omit(unstandardised$stats$lat_r))), type = "l", lty = 2, col = 2)
  axis(1, labels = F)
  axis(4)
  mtext(side = 4, line = 3, "Latitude range (deg)", cex = 0.8)
  
  par(mar = c(4.5, 4.5, 0.5, 0.5))
  
  plot(mpt, unlist(lapply(unstandardised$data, nrow)), xlim = c(stg[1], stg[length(stg)]), ylim = c(0, max(unlist(lapply(unstandardised$data, nrow)))), xlab = "Time (Ma)", type = "l", ylab = "Occurrences (n)")
  lines(mpt, unlist(lapply(standardised$data, nrow)), lty = 2)
  
  par(mar = c(4.5, 4.5, 0.5, 4.5))
  
  if(all(is.na(unstandardised$stats$mst_l))) {
    plot(NA, NA, xlim = c(stg[1], stg[length(stg)]), xlab = "Time (Ma)", ylab = "MST length (km)", ylim = c(0, 0))
  } else {
    plot(mpt, standardised$stats$mst_l, xlim = c(stg[1], stg[length(stg)]), xlab = "Time (Ma)", ylab = "MST length (km)", type = "l", lty = 2)
  }
  par(new = TRUE)
  plot(mpt, standardised$div$qD_0.5, xlim = c(stg[1], stg[length(stg)]), axes = FALSE, type = "l", col = 3, ylab = "", xlab = "")
  axis(4)
  mtext(side = 4, line = 3, "SQS diversity (q = 0.5)", cex = 0.8)
  dev.off()
  
  pv <- c(cor.test(standardised$div[,2], standardised$stats$mst_l, method = "pearson", alternative = "greater")$p.value,
          cor.test(standardised$div[,2], standardised$stats$mst_l, method = "spearman", alternative = "greater")$p.value,
          cor.test(standardised$div[,2], standardised$stats$lng_r, method = "pearson", alternative = "greater")$p.value,
          cor.test(standardised$div[,2], standardised$stats$lng_r, method = "spearman", alternative = "greater")$p.value,
          cor.test(standardised$div[,2], standardised$stats$lat_r, method = "pearson", alternative = "greater")$p.value,
          cor.test(standardised$div[,2], standardised$stats$lat_r, method = "spearman", alternative = "greater")$p.value)
  
  rv <- c(cor.test(standardised$div[,2], standardised$stats$mst_l, method = "pearson", alternative = "greater")$estimate,
          cor.test(standardised$div[,2], standardised$stats$mst_l, method = "spearman", alternative = "greater")$estimate,
          cor.test(standardised$div[,2], standardised$stats$lng_r, method = "pearson", alternative = "greater")$estimate,
          cor.test(standardised$div[,2], standardised$stats$lng_r, method = "spearman", alternative = "greater")$estimate,
          cor.test(standardised$div[,2], standardised$stats$lat_r, method = "pearson", alternative = "greater")$estimate,
          cor.test(standardised$div[,2], standardised$stats$lat_r, method = "spearman", alternative = "greater")$estimate)
  
  write.table(cbind.data.frame(rv = round(rv, digits = 3), pv = round(pv, digits = 3)), paste0("tables/", nms[s], "_corr.txt"),
              row.names = c("mst_l_pearson", "mst_l_spearman", "lng_r_pearson", "lng_r_spearman", "lat_r_pearson", "lat_r_spearman"), quote = F)
  saveRDS(standardised, paste0("data/", nms[s], ".RData"))
  
  standardised_unstd <- spacetimestand(x = occs, lng = "Paleo_Longitude", lat = "Paleo_Latitude", area = alist[[s]],
                                 mst = TRUE, mst_l = amst_unstd[s], lnglat = alnglat[s], lng_r = lng_unstd[s], lat_r = lat_unstd[s],
                                 intervals = stg, max_ma = "Max_Ma", min_ma = "Min_Ma", div = T, tax = "Taxon_name")
  
  pv <- c(cor.test(standardised_unstd$div[,2], standardised_unstd$stats$mst_l, method = "pearson", alternative = "greater")$p.value,
          cor.test(standardised_unstd$div[,2], standardised_unstd$stats$mst_l, method = "spearman", alternative = "greater")$p.value,
          cor.test(standardised_unstd$div[,2], standardised_unstd$stats$lng_r, method = "pearson", alternative = "greater")$p.value,
          cor.test(standardised_unstd$div[,2], standardised_unstd$stats$lng_r, method = "spearman", alternative = "greater")$p.value,
          cor.test(standardised_unstd$div[,2], standardised_unstd$stats$lat_r, method = "pearson", alternative = "greater")$p.value,
          cor.test(standardised_unstd$div[,2], standardised_unstd$stats$lat_r, method = "spearman", alternative = "greater")$p.value)
  
  rv <- c(cor.test(standardised_unstd$div[,2], standardised_unstd$stats$mst_l, method = "pearson", alternative = "greater")$estimate,
          cor.test(standardised_unstd$div[,2], standardised_unstd$stats$mst_l, method = "spearman", alternative = "greater")$estimate,
          cor.test(standardised_unstd$div[,2], standardised_unstd$stats$lng_r, method = "pearson", alternative = "greater")$estimate,
          cor.test(standardised_unstd$div[,2], standardised_unstd$stats$lng_r, method = "spearman", alternative = "greater")$estimate,
          cor.test(standardised_unstd$div[,2], standardised_unstd$stats$lat_r, method = "pearson", alternative = "greater")$estimate,
          cor.test(standardised_unstd$div[,2], standardised_unstd$stats$lat_r, method = "spearman", alternative = "greater")$estimate)
  
  write.table(cbind.data.frame(rv = round(rv, digits = 3), pv = round(pv, digits = 3)), paste0("tables/", nms[s], "_corr_unstd.txt"),
              row.names = c("mst_l_pearson", "mst_l_spearman", "lng_r_pearson", "lng_r_spearman", "lat_r_pearson", "lat_r_spearman"), quote = F)
}


Global <- readRDS("data/Global.RData")
NAfEu <- readRDS("data/NAfEu.RData")
WAfSAm <- readRDS("data/WAfSAm.RData")
stgn <- c("Bajocian–Callovian", "Oxfordian", "Kimmeridgian", "Tithonian", "Berriasian", "Valanginian",
          "Hauterivian", "Barremian", "Aptian", "Albian", "Cenomanian–Turonian")
roccs_s <- lapply(list(Global, NAfEu, WAfSAm), function(x) {unique(unlist(lapply(x$data, function(y) {y$Occurrence_Number})))})
names(roccs_s) <- c("Global", "NAfEu", "WAfSAm")

occs$globe <- occs$region <- NA
occs$globe[which(occs$Occurrence_Number %in% unlist(roccs_s$Global))] <- "Global"
occs$region[which(occs$Occurrence_Number %in% unlist(roccs_s$NAfEu))] <- "NAfEu"
occs$region[which(occs$Occurrence_Number %in% unlist(roccs_s$WAfSAm))] <- "WAfSAm"

tiff("map.tiff", width = 2048, height = 1250, res = 300)

par(mfrow = c(3, 3), mar = c(0, 0, 0, 0))  

for (i in 1:length(shelf_list)) {
  binned <- occs[which(occs$Max_Ma <= stg[i] & occs$Max_Ma >= stg[i + 1] |
                         occs$Min_Ma <= stg[i] & occs$Min_Ma >= stg[i + 1]), ]
  
  plot(frame, col = "aliceblue", axes = FALSE, xlab = "", ylab = "", lwd = 0.7)
  plot(shelf_list[[i]], col = "lightsteelblue1", border = "lightsteelblue2", add = TRUE, lwd = 0.7)
  plot(coast_list[[i]], col = "seashell", border = "steelblue2", add = TRUE, lwd = 0.7)
  plot(frame, add = TRUE)
  
  point_colors <- ifelse(binned$Original_Assignment %in% c("WAf", "SAm"),"#458B00",
                                ifelse(binned$Original_Assignment %in% c("NAf", "Eu"), "purple",
                                              "#CDC8B1"))
  
  points(binned$Paleo_Longitude, binned$Paleo_Latitude,
         pch = 4, 
         bg = adjustcolor(point_colors),
         col = point_colors,
         cex = 0.5)
}

dev.off()

write.table(occs, "data/Cyprideidae_occ_standardized.txt", sep = "\t", row.names = FALSE, quote = FALSE)

