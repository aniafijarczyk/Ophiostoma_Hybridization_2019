library(ggplot2)
library(ggrepel)
library(dplyr)
library(grid)



# Function: get the x and y positions of a single ggrepel label
get.xy.pos.labs <- function(n) {
  grb <- grid.get(n)
  data.frame(
    x = xrg[1]+diff(xrg)*convertX(grb$x, "native", valueOnly = TRUE),
    y = yrg[1]+diff(yrg)*convertY(grb$y, "native", valueOnly = TRUE)
  )
}


# Reading dataset
d <- read.table("samples_locations.txt", sep = "\t", header = TRUE, na.strings = "NA")
head(d)
x = d$Longitude
y = d$Latitude
ShortSci = d$ID
df <- data.frame(x = x, y = y, z = as.factor(ShortSci))




# Running simulations for N America and Euroasia

fd <- df %>% filter(z != '81' & z != '82' & z != '28')

p1 <- ggplot(data = fd, aes(x = x, y = y)) + theme_bw() + 
  scale_x_continuous(limits = c(-150, 110)) +
  scale_y_continuous(limits = c(30, 61)) +
  geom_text_repel(aes(label = z), 
                  box.padding = unit(0.1, "lines"),
                  nudge_y = 3,
                  nudge_x = 1,
                  force=0.9,
                  point.padding = 0.1) +
  geom_point(colour = "green", size = 3)
p1

## Get x and y plot ranges
ggp1 <- ggplot_build(p1)
xrg <- ggp1$layout$panel_params[[1]]$x.range
yrg <- ggp1$layout$panel_params[[1]]$y.range

grid.force()
kids <- childNames(grid.get("textrepeltree", grep = TRUE))

## Get positions of all ggrepel labels
dts <- do.call(rbind, lapply(kids, get.xy.pos.labs))
dts$lab <- fd$z
print(dts)



# Running simulations for New Zealand

fz <- df %>% filter(z == '81' | z == '82')

p2 <- ggplot(data = fz, aes(x = x, y = y)) + theme_bw() + 
  scale_x_continuous(limits = c(150, 200)) +
  scale_y_continuous(limits = c(-50, -20)) +
  geom_text_repel(aes(label = z), 
                  box.padding = unit(0.6, "lines"),
                  nudge_y = 3) +
  geom_point(colour = "green", size = 3)
p2

## Get x and y plot ranges
ggp2 <- ggplot_build(p2)
xrg <- ggp2$layout$panel_params[[1]]$x.range
yrg <- ggp2$layout$panel_params[[1]]$y.range

grid.force()
kids2 <- childNames(grid.get("textrepeltree", grep = TRUE))

## Get positions of all ggrepel labels
dts2 <- do.call(rbind, lapply(kids2, get.xy.pos.labs))
dts2$lab <- fz$z
print(dts2)



# Concatenating two datasets

DTS <- rbind(dts, dts2)

write.table(DTS,file = "simulated_coords.txt", sep="\t", col.names = TRUE, row.names = FALSE)
