setwd("/media/anna/Volume/Ophiostoma/pop_genomics/2019_08_paper_scripts/Strcuture_and_map_O_novo_ulmi")
source('Pipie.R', chdir = TRUE)
library(dplyr)
library(tidyr)
library(plotrix)
library(maps)
library(mapdata)
library(maptools)

### Importing file with borders
fdc <- read.table("world_intermediate.txt", header = F)

### Importing file containing the coordinate of strains
df <- read.table("samples_locations.txt", sep = "\t", header = TRUE, na.strings = "NA")
s <- read.table("simulated_coords.txt", sep = "\t", header = TRUE, col.names = c("Long_sim", "Lat_sim", "lab"))
m <- merge(s, df, by.x = "lab", by.y = "ID", sort = FALSE)

### Importing file containing the output of structure
K3 <- read.table("k3.csv", sep = ",")
head(K3)

### Setting colors
cl <- colors()
col_ame1 <- "#d55e00"
col_ame2 <- "#f0e442"
col_nov <- "#009e73"
col_ou <- "#0072b2"

### Plotting the map of Northern Hemispere


d <- m %>% filter(lab != '81' & lab != '82')

pdf("Figure_piechart_novoulmi.pdf", 12, 6.5)

map("world",xlim = c(-140, 110), ylim = c(30, 70), fill = F, col = colors()[200], type = "l")
map.axes(xaxt = "n", yaxt = "n")
for(i in as.character(d$Isolate)){
  di <- filter(d, Isolate == i)
  K3i <- filter(K3, V2 == i) %>% dplyr::select(-V1, -V2) %>% as.matrix
  if (length(K3i) > 0) {
    points(di$Longitude, di$Latitude, pch = 21, col = "black", bg = "black", cex = 0.4)
    lines(x = c(di$Long_sim, di$Longitude), y = c(di$Lat_sim, di$Latitude), col="black")
  }
}
for(i in as.character(d$Isolate)){
	# Extract the strain to obtain its location
	di <- filter(d, Isolate == i)
	# Extract the strain to obtain its results in structure
	K3i <- filter(K3, V2 == i) %>% dplyr::select(-V1, -V2) %>% as.matrix
	# If the strain is indeed in the stroucture output ...
	if (length(K3i) > 0) {
		pipie(K3i, long = di$Long_sim, lat = di$Lat_sim, labels = "", radius = 0.05, 
		      col = c(col_ame1, col_ame2, col_nov),border=FALSE)
	}
}
dev.off()


### Plotting the New Zealand

d <- m %>% filter(lab == '81' | lab == '82')

pdf("Figure_piechart_novoulmi_NZ.pdf", 2.6, 5.2)

{
  map("world",xlim = c(166, 181), ylim = c(-48, -31), fill = F, col = colors()[200], type = "l", bg = "white")
  map.axes(xaxt = "n", yaxt = "n")
  points(d$Longitude, d$Latitude, pch = 21, col = "black", bg = "black", cex = 0.4)

  # For each strains ...
  for(i in as.character(d$Isolate)){
  	# Extract the strain to obtain its location	
  	di <- filter(d, Isolate == i)
  	# Extract the strain to obtain its results in structure
  	K3i <- filter(K3, V2 == i) %>% select(-V1, -V2) %>% as.matrix
  
  	# If the strain is indeed in the stroucture output ...
  	if (length(K3i) > 0) {
  		# print a psmall black point at the sampling location
  		points(di$Longitude, di$Latitude, pch = 21, col = "black", bg = "black", cex = 0.4)
  		# print a line that will link the sampling location with the pie chart (structure results of the strain)
  		lines(x = c(di$Long_sim, di$Longitude), y = c(di$Lat_sim, di$Latitude), col="black")
  		# print the pie chart	
  		pipie(K3i, long = di$Long_sim, lat = di$Lat_sim, labels = "",
  		      radius = 0.09, col = c(col_ame1, col_ame2, col_nov),border=FALSE)
  	}
  }
}

dev.off()





