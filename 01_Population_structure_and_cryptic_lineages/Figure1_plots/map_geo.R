

library(plotrix) # pour les labels
library(dplyr)
library(maps)
library(mapdata)


# Setting colors

novo_ulmi = "olivedrab4"
americana = "orange3"
ulmi = "dodgerblue3"
himal_ulmi = "mediumpurple4" 
ame1 = "firebrick3"
ame2 = "goldenrod2"

# Reading maps and samples

fdc <- read.table("world_intermediate.txt", header = F)
d <- read.table("Liste_finale_plaque1_LN2.txt", header = T)



# WORLD

pdf("World_North_col_AmeDistinction_blackWhite.pdf", 17, 6.5)
{
  map("world",xlim = c(-140, 100), ylim = c(30, 70), fill = F, col = colors()[200], type = "l")
  map.axes(xaxt = "n", yaxt = "n")
  
  points(d$Longitude, d$Latitude, pch = 21, col = "black", bg = "black", cex = 0.4)
  
  for(i in 1:nrow(d)){
    lines(x = c(d$Long_sim[i], d$Longitude[i]), y = c(d$Lat_sim[i], d$Latitude[i]), col="black")
    if (d[i,30] == "NOV") {points( x = d[i,9], y = d[i,8], pch = 21, cex=2, col=novo_ulmi, bg = novo_ulmi)}
    if (d[i,30] == "AME1") {points( x = d[i,9], y = d[i,8], pch = 22, cex=2, col=ame1, bg = ame1)}
    if (d[i,30] == "AME2") {points( x = d[i,9], y = d[i,8], pch = 22, cex=2, col=ame2, bg = ame2)}
    if (d[i,30] == "ULM") {points( x = d[i,9], y = d[i,8], pch = 24, cex=2, col=ulmi, bg = ulmi)}
    if (d[i,30] == "HIM") {points( x = d[i,9], y = d[i,8], pch = 23, cex=2, col=himal_ulmi, bg = himal_ulmi)}
  }

  legend("topleft", legend = c("AME1", "AME2", "O. himal-ulmi","NOV", "ULM"), pch = c(22, 22, 23, 21, 24), pt.bg = c(ame1, ame2, himal_ulmi, novo_ulmi, ulmi),pt.cex=2, cex = 1, ncol = 2, bg = "white", text.font = c(1,1,3,1,1))
}
dev.off()





### New Zaeland

pdf("World_NZ_col_AmeDistinction_BlackWhite_noborders.pdf", 2.2, 5.2)
{
  #plot(fdc, type = "l", col = "gray60", xlab = "long", ylab = "lat", xlim = c(165, 180), ylim = c(-50, -25))
  map("world",xlim = c(166, 181), ylim = c(-48, -31), fill = F, col = colors()[200], type = "l", bg = "white")
  #map.axes()
  map.axes(xaxt = "n", yaxt = "n")
  points(d$Longitude, d$Latitude, pch = 21, col = "black", bg = "black", cex = 0.4)
  for(i in 1:nrow(d)){
    lines(x = c(d$Long_sim[i], d$Longitude[i]), y = c(d$Lat_sim[i], d$Latitude[i]), col="black")
    if (d[i,30] == "NOV") {points( x = d[i,9], y = d[i,8], pch = 21, cex=2, col=novo_ulmi, bg = novo_ulmi)}
    if (d[i,30] == "AME1") {points( x = d[i,9], y = d[i,8], pch = 22, cex=1.75, col=ame1, bg = ame1)}
    if (d[i,30] == "AME2") {points( x = d[i,9], y = d[i,8], pch = 22, cex=2, col=ame2, bg = ame2)}
    if (d[i,30] == "ULM") {points( x = d[i,9], y = d[i,8], pch = 24, cex=2, col=ulmi, bg = ulmi)}
    if (d[i,30] == "HIM") {points( x = d[i,9], y = d[i,8], pch = 23, cex=2, col=himal_ulmi, bg = himal_ulmi)}
  }
}
dev.off()

