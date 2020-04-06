################### World map

library(rworldmap)
library(dplyr)
library(ggplot2)
library(geosphere)
library(gpclib)



worldMap <- getMap()
world.points <- fortify(worldMap)
world.points$region <- world.points$id

world.df <- world.points[,c("long","lat","group", "region")]

worldmap <- ggplot() + 
  geom_polygon(data = world.df, aes(x = long, y = lat, group = group)) +
  scale_y_continuous(breaks = (-2:2) * 30) +
  scale_x_continuous(breaks = (-4:4) * 45)

worldmap

pdf("plotGlobe.pdf",w=8,h=8,compress=FALSE)

worldmap <- ggplot() + 
  geom_polygon(data = world.df, aes(x = long, y = lat, group = group)) +
  scale_y_continuous(breaks = (-2:2) * 30) +
  #scale_x_continuous(breaks = (-4:4) * 45) +
  scale_x_continuous(breaks = c(55,100,145,190,-135)) +
  scale_fill_manual(values = "red") +
  coord_map("ortho", orientation=c(0, 140, 0)) +
  theme_bw() +
  theme(panel.border = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank())
worldmap


dev.off()
