# 17 april 2018
# Strucutre pie chart on map
setwd("/media/anna/Volume/Ophiostoma/pop_genomics/2019_04_paper_figures/Structure_pie_map")
library(dplyr)
source('Pipie.R', chdir = TRUE)


# File containing the map
fdc <- read.table("world_intermediate.txt", header = F)

# File containing the coordinate of strains
d <- read.table("Liste_finale_plaque1_LN.txt", header = T)
head(d)
str(d)

# File containing the output of structure
K5 <- read.table("k5.csv", sep = ",")
head(K5)

# a small function to handle more easily the color
cl <- colors()

################################
# Plot of the North Hemisphere #
################################

# Plot the map
plot(fdc, type = "l", col = "gray60", xlab = "long", ylab = "lat", xlim = c(-140, 100), ylim = c(30, 70))


# For each strains ...
for(i in as.character(d$Isolate)){
	#i="DDS100"
	
	# Extract the strain to obtain its location
	di <- filter(d, Isolate == i)
	
	# Extract the strain to obtain its results in structure
	K5i <- filter(K5, V2 == i) %>% select(-V1, -V2) %>% as.matrix

	# If the strain is indeed in the stroucture output ...
	if (length(K5i) > 0) {
		
		# print a psmall black point at the sampling location
		points(di$Longitude, di$Latitude, pch = 21, col = "black", bg = "black", cex = 0.4)
		
		# print a line that will link the sampling location with the pie chart (structure results of the strain)
		lines(x = c(di$Long_sim, di$Longitude), y = c(di$Lat_sim, di$Latitude), col="black")
		
		# print the pie chart	
		pipie(K5i, long = di$Long_sim, lat = di$Lat_sim, labels = "", radius = 0.05, col = c(cl[517], cl[149], cl[497], cl[136], cl[12]))
	}
	
}



###########################
# Plot of the New Zealand #
###########################

# Plot the map
plot(fdc, type = "l", col = "gray60", xlab = "long", ylab = "lat", xlim = c(165, 180), ylim = c(-50, -25))


# For each strains ...
for(i in as.character(d$Isolate)){
	#i="DDS100"
	
	# Extract the strain to obtain its location	
	di <- filter(d, Isolate == i)

	# Extract the strain to obtain its results in structure
	K5i <- filter(K5, V2 == i) %>% select(-V1, -V2) %>% as.matrix

	# If the strain is indeed in the stroucture output ...
	if (length(K5i) > 0) {

		# print a psmall black point at the sampling location
		points(di$Longitude, di$Latitude, pch = 21, col = "black", bg = "black", cex = 0.4)
		
		# print a line that will link the sampling location with the pie chart (structure results of the strain)
		lines(x = c(di$Long_sim, di$Longitude), y = c(di$Lat_sim, di$Latitude), col="black")
	
		# print the pie chart	
		pipie(K5i, long = di$Long_sim, lat = di$Lat_sim, labels = "", radius = 0.09, col = c(cl[517], cl[149], cl[497], cl[136], cl[12]))
	}
	
}






