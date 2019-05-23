
# Script for plotting color matrices for estimates of divergence and diversity in Ophiostoma lineages

library(ggplot2)
library(gridExtra)

### Plot of global pairwise distance (dxy) and diversity

d <- read.table(file = "matrix_global_estimates.txt",header=TRUE,sep="\t")


p1 <- ggplot(d) + 
  aes(x = SP2, y = SP1, fill = Value) +
  geom_tile(color = "white") +
  geom_text(aes(SP2, SP1, label = format(Value, scientific = FALSE)), color = "black", size = 4) +
  scale_fill_gradient(low = "white", high = "darkgrey", 
                      limit = c(0,0.04), space = "Lab",
                      name="Global\npairwise\ndistance") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(vjust = 1, size = 12, hjust = 1),
        axis.text.y = element_text(size = 12),
        axis.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        #legend.justification = c(1, 0),
        legend.position = c(0.25, 0.8),
        legend.direction = "horizontal") +
  guides(fill = guide_colorbar(barwidth = 8, barheight = 1.5,
                               title.position = "top", title.hjust = 0.5)) +
  coord_fixed()
p1



### Plot of per-window pairwise distance (dxy) and diversity

w <- read.table(file = "matrix_window_estimates.txt",header=TRUE,sep="\t")

p2 <- ggplot(w) + 
  aes(x = SP2, y = SP1, fill = Value) +
  geom_tile(color = "white") +
  geom_text(aes(SP2, SP1, label = format(Value, scientific = FALSE)), color = "black", size = 4) +
  scale_fill_gradient(low = "white", high = "darkgrey", 
                      limit = c(0,0.04), space = "Lab",
                      name="Per window\npairwise\ndistance") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(vjust = 1, size = 12, hjust = 1),
        axis.text.y = element_text(size = 12),
        axis.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        #legend.justification = c(1, 0),
        legend.position = c(0.25, 0.8),
        legend.direction = "horizontal") +
  guides(fill = guide_colorbar(barwidth = 8, barheight = 1.5,
                               title.position = "top", title.hjust = 0.5)) +
  coord_fixed()
p2


pdf("plotMatrix.pdf", w=10, h=4)
grid.arrange(p1,p2,ncol=2)
dev.off()
