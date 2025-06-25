library(palaeoverse)
library(dplyr)

paleocoord <- function(input_file, output_file) {
  modern_coord <- read.table(input_file, header = TRUE) %>%
    mutate(Midpoint_age = (Min_Ma + Max_Ma)/2) %>%
    select(Longitude, Latitude, Midpoint_age)
  res_data <- palaeorotate(
    modern_coord,
    lng = "Longitude",
    lat = "Latitude",
    age = "Midpoint_age",
    model = "PALEOMAP"
  ) %>% select(p_lat, p_lng)
  write.table(res_data, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)
}


paleocoord("North_Africa_coord.txt", "North_Africa_coord_PALEOMAP.txt")
paleocoord("West_Africa_coord.txt", "West_Africa_coord_PALEOMAP.txt")
paleocoord("North_America_coord.txt", "North_America_coord_PALEOMAP.txt")
paleocoord("Europe_coord.txt", "Europe_coord_PALEOMAP.txt")
paleocoord("South_America_coord.txt", "South_America_coord_PALEOMAP.txt")
