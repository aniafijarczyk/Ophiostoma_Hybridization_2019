library(plotrix)
library(dplyr)
library(maps)
library(mapdata)
library(ggplot2)
library(maptools)

# Setting colors

novo_ulmi = "olivedrab4"
americana = "orange3"
ulmi = "dodgerblue3"
himal_ulmi = "mediumpurple4" 
ame1 = "firebrick3"
ame2 = "goldenrod2"

# Reading maps and samples

fdc <- read.table("world_intermediate.txt", header = T, col.names = c("Lat","Long"))
head(fdc)

# Reading cordinates and simulated coordinates

df <- read.table("samples_locations.txt", sep = "\t", header = TRUE, na.strings = "NA")
s <- read.table("simulated_coords.txt", sep = "\t", header = TRUE, col.names = c("Long_sim", "Lat_sim", "lab"))
m <- merge(s, df, by.x = "lab", by.y = "ID", sort = FALSE)

### Plotting North America and Eurasia

d <- m %>% filter(lab != '81' & lab != '82')

pdf("World_North_col.pdf", 17, 6.5)
{
  map("world",xlim = c(-140, 110), ylim = c(30, 70), fill = F, col = colors()[200], type = "l")
  map.axes(xaxt = "n", yaxt = "n")
  points(d$Longitude, d$Latitude, pch = 21, col = "black", bg = "black", cex = 0.4)
 
  for(i in 1:nrow(d)){
    lines(x = c(d$Long_sim[i], d$Longitude[i]), y = c(d$Lat_sim[i], d$Latitude[i]), col="black")
  }
  for(i in 1:nrow(d)){
    if (d$Lineage[i] == "NOV") {points( x =  d$Long_sim[i], y = d$Lat_sim[i], pch = 21, cex=1.7, col=novo_ulmi, bg = novo_ulmi)}
    if (d$Lineage[i] == "AME1") {points( x = d$Long_sim[i], y = d$Lat_sim[i], pch = 22, cex=1.7, col=ame1, bg = ame1)}
    if (d$Lineage[i] == "AME2") {points( x = d$Long_sim[i], y = d$Lat_sim[i], pch = 22, cex=1.7, col=ame2, bg = ame2)}
    if (d$Lineage[i] == "OU") {points( x = d$Long_sim[i], y = d$Lat_sim[i], pch = 24, cex=1.7, col=ulmi, bg = ulmi)}
    if (d$Lineage[i] == "OHU") {points( x = d$Long_sim[i], y = d$Lat_sim[i], pch = 23, cex=1.7, col=himal_ulmi, bg = himal_ulmi)}
  }
  
  legend("topleft", legend = c("AME1", "AME2", "NOV", "OU", "O. himal-ulmi"), pch = c(22, 22, 21, 24, 23), 
         pt.bg = c(ame1, ame2, novo_ulmi, ulmi, himal_ulmi),pt.cex=2, cex = 1.4, ncol = 5, bg = "white", text.font = c(1,1,1,1,3))
}
dev.off()





### Plotting New Zaeland

d <- f %>% filter(lab == '81' | lab == '82')

pdf("World_NZ_col.pdf", 2.2, 5.2)
{
  map("world",xlim = c(166, 181), ylim = c(-48, -31), fill = F, col = colors()[200], type = "l", bg = "white")
  map.axes(xaxt = "n", yaxt = "n")
  points(d$Longitude, d$Latitude, pch = 21, col = "black", bg = "black", cex = 0.4)
  
  for(i in 1:nrow(d)){
    lines(x = c(d$Long_sim[i], d$Longitude[i]), y = c(d$Lat_sim[i], d$Latitude[i]), col="black")
  }
  for(i in 1:nrow(d)){
    if (d$Lineage[i] == "NOV") {points( x =  d$Long_sim[i], y = d$Lat_sim[i], pch = 21, cex=1.7, col=novo_ulmi, bg = novo_ulmi)}
    if (d$Lineage[i] == "AME1") {points( x = d$Long_sim[i], y = d$Lat_sim[i], pch = 22, cex=1.7, col=ame1, bg = ame1)}
    if (d$Lineage[i] == "AME2") {points( x = d$Long_sim[i], y = d$Lat_sim[i], pch = 22, cex=1.7, col=ame2, bg = ame2)}
    if (d$Lineage[i] == "OU") {points( x = d$Long_sim[i], y = d$Lat_sim[i], pch = 24, cex=1.7, col=ulmi, bg = ulmi)}
    if (d$Lineage[i] == "OHU") {points( x = d$Long_sim[i], y = d$Lat_sim[i], pch = 23, cex=1.7, col=himal_ulmi, bg = himal_ulmi)}
  }
}
dev.off()

