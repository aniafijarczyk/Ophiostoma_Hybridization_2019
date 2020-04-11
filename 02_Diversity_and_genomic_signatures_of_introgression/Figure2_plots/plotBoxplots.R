library(ggplot2)
library(dplyr)
library(gridExtra)
library(grid)


# Reading data

df1 <- read.table("./data/ame1.ame2_circos_Dxy.txt", sep = "\t", header = FALSE, col.names = c("chrom","start","stop","stat"))
df2 <- read.table("./data/nov.ame1_circos_Dxy.txt", sep = "\t", header = FALSE, col.names = c("chrom","start","stop","stat"))
df3 <- read.table("./data/nov.ame2_circos_Dxy.txt", sep = "\t", header = FALSE, col.names = c("chrom","start","stop","stat"))
df4 <- read.table("./data/ulm.ame1_circos_Dxy.txt", sep = "\t", header = FALSE, col.names = c("chrom","start","stop","stat"))
df5 <- read.table("./data/ulm.ame2_circos_Dxy.txt", sep = "\t", header = FALSE, col.names = c("chrom","start","stop","stat"))
df6 <- read.table("./data/ulm.nov_circos_Dxy.txt", sep = "\t", header = FALSE, col.names = c("chrom","start","stop","stat"))

df7 <- read.table("./data/ame1_circos_tP.txt", sep = "\t", header = FALSE, col.names = c("chrom","start","stop","stat"))
df8 <- read.table("./data/ame2_circos_tP.txt", sep = "\t", header = FALSE, col.names = c("chrom","start","stop","stat"))
df9 <- read.table("./data/nov_circos_tP.txt", sep = "\t", header = FALSE, col.names = c("chrom","start","stop","stat"))
df10 <- read.table("./data/ulm_circos_tP.txt", sep = "\t", header = FALSE, col.names = c("chrom","start","stop","stat"))


# Preparing data frames

df7["sp"] <- "AME1"
df8["sp"] <- "AME2"
df9["sp"] <- "NOV"
df10["sp"] <- "OU"

df1["sp"] <- "AME1-AME2"
df2["sp"] <- "AME1-NOV"
df3["sp"] <- "AME2-NOV"

df4["sp"] <- "AME1-OU"
df5["sp"] <- "AME2-OU"
df6["sp"] <- "NOV-OU"

df_pi <- rbind(df7,df8,df9,df10)
dxy_intra <- rbind(df1,df2,df3)
dxy_inter <- rbind(df4,df5,df6)


# Creating dataframes with mean values

pi <- data.frame("sp" = c("AME1","AME2","NOV","OU"),
                 "stat" = c(0.0016,0.0004,0.0018,0.0028),
                 "sd" = c(0.0007,0.0002,0.0008,0.0004))

dxy1 <- data.frame("sp" = c("AME1-AME2","AME1-NOV","AME2-NOV"),
                   "stat" = c(0.0027,0.0051,0.0053))

dxy2 <- data.frame("sp" = c("AME1-OU","AME2-OU", "NOV-OU"),
                   "stat" = c(0.0323,0.0324,0.0319))

pal1 <- wes_palette("Royal2")
col1 <- pal1[1]
col1 <- "grey90"


# Plotting boxplots

p1 <- ggplot(df_pi) + aes(x = sp, y = stat, fill = sp) +
  scale_y_continuous(limits = c(0,0.05)) +
  labs(x = "Nucleotide diversity", y = "pairwise distance") +
  scale_fill_manual(values = c("#cd2626","#eeb422","#698b22","#1874cd","#5d478b")) +
  geom_boxplot() +
  #geom_rect(xmin = 0.6,xmax = 1.4, ymin = 0.0016, ymax = 0.0016, color = "red") +
  geom_point(data = pi, aes(x = sp, y = stat), shape = 23, fill = col1, colour = "black", size = 3) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.text.x = element_text(size = 15, angle = 45, hjust =1),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_blank(),
        legend.position = 'None')
p1

p2 <- ggplot(dxy_intra) + aes(x = sp, y = stat, fill = sp) +
  scale_y_continuous(limits = c(0,0.05)) +
  geom_boxplot() +
  scale_fill_manual(values = c("grey90","grey90","grey90")) +
  geom_point(data = dxy1, aes(x = sp, y = stat), shape = 23, fill = col1, colour = "black", size = 3) +
  theme(axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.text.x = element_text(size = 15, angle = 45, hjust =1),
        axis.title.x = element_blank(),
        legend.position = 'None')
p2

p3 <- ggplot(dxy_inter) + aes(x = sp, y = stat, fill =sp) +
  scale_y_continuous(limits = c(0,0.05)) +
  geom_boxplot() +
  scale_fill_manual(values = c("grey70","grey70","grey70")) +
  labs(x = "Divergence between") +
  geom_point(data = dxy2, aes(x = sp, y = stat), shape = 23, fill = col1, colour = "black", size = 3) +
  theme(axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.text.x = element_text(size = 15, angle = 45, hjust =1),
        legend.position = 'None')
p3


# Combining plots

gt1 <- ggplotGrob(p1)
gt2 <- ggplotGrob(p2)
gt3 <- ggplotGrob(p3)
newWidth = unit.pmax(gt1$widths[2:3], gt2$widths[2:3],gt3$widths[2:3])
gt1$widths[2:3] = as.list(newWidth)
gt2$widths[2:3] = as.list(newWidth)
gt3$widths[2:3] = as.list(newWidth)

newHeight = unit.pmax(gt1$heights[8:9], gt2$heights[8:9],gt3$heights[8:9])
gt1$heights[8:9] = as.list(newHeight)
gt2$heights[8:9] = as.list(newHeight)
gt3$heights[8:9] = as.list(newHeight)


png("plotBoxplots.png",w=7,h=5,res=150,units="in")
grid.arrange(gt1, gt2, gt3, ncol=3)
dev.off()

