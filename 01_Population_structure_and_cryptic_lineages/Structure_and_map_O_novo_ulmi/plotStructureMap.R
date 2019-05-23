
library(dplyr)
source('Pipie.R', chdir = TRUE)


### Importing file containing the map
fdc <- read.table("world_intermediate.txt", header = F)

### Importing file containing the coordinate of strains
d <- read.table("Liste_finale_plaque1_LN2.txt", header = T)
head(d)
str(d)

### Importing file containing the output of structure
K3 <- read.table("k3.csv", sep = ",")
head(K3)

### A small function to handle more easily the color
cl <- colors()
goldenrod2 <- colors()[149]
olivedrab4 <- colors()[497]
firebrick3 <- colors()[136]


### Plotting the map of Northern Hemispere


#pdf("Figure_piechart_novoulmi.pdf", 8.82, 4.60)

plot(fdc, type = "l", col = "gray60", xlab = "Longitude", ylab = "Latitude", xlim = c(-140, 100), ylim = c(30, 70))

# For each strain ...
for(i in as.character(d$Isolate)){
	# i="DDS100"
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
		pipie(K3i, long = di$Long_sim, lat = di$Lat_sim, labels = "", radius = 0.05, 
		      col = c(cl[136], cl[149], cl[497]),border=FALSE)
	}
}

#dev.off()


### Plotting the New Zealand

#pdf("Figure_piechart_novoulmi_NZ.pdf", 2.6, 5.2)

# Plot the map
plot(fdc, type = "l", col = "gray60", xlab = "", ylab = "", xlim = c(165, 180), ylim = c(-50, -25))

# For each strains ...
for(i in as.character(d$Isolate)){
	#i="DDS100"
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
		      radius = 0.09, col = c(cl[136], cl[149], cl[497]),border=FALSE)
	}
}

#dev.off()





