library(fields)
library(maps)
# Adjusting genotype-environment association models
library(LEA)
library(tidyverse)
library(vcfR)
library(SNPfiltR)
library(LEA)

genotype = LEA::read.geno("/users/l/e/lericsso/projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.thinned.geno")
coordinates = as.matrix(read.table(meta2[,7:8]))
latitude <- meta2$latitude
longitude <- meta2$longitude

plot(longitude, latitude, cex = .4, col = "darkblue",
     xlab = "Longitude", ylab = "Latitude",
     main = "coordinates", las = 1)
maps::map(add = TRUE, interior = FALSE, col = "grey40")

########################


# loading the required packages
library(ggplot2)
library(ggmap)

# creating a sample data.frame with your lat/lon points
lon <- c(meta2$longitude)
lat <- c(meta2$latitude)
df <- as.data.frame(cbind(lon,lat))

# getting the map
mapgilbert <- get_stadiamap(location = c(lon = mean(df$lon), lat = mean(df$lat)), zoom = 4,
                      maptype = "satellite", scale = 2)

# plotting the map with some points on it
ggmap(mapgilbert) +
  geom_point(data = df, aes(x = lon, y = lat, fill = "red", alpha = 0.8), size = 5, shape = 21) +
  guides(fill=FALSE, alpha=FALSE, size=FALSE)


