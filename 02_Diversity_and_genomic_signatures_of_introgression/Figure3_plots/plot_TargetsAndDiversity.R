
library(ggplot2)
library(gridExtra)


### Reading files and plotting diversity

ame1 <- read.table("combineTables_AME1.txt", sep = "\t", header = TRUE)
ame2 <- read.table("combineTables_AME2.txt", sep = "\t", header = TRUE)
nov <- read.table("combineTables_NOV.txt", sep = "\t", header = TRUE)

a1 <- ame1[c("Chr","WinCenter","Pi_AME1","Class","Win","CommonAME1")]
a2 <- ame2[c("Chr","WinCenter","Pi_AME2","Class","Win","CommonAME2")]
n <- nov[c("Chr","WinCenter","Pi_NOV","Class","Win","CommonNOV")]

colnames(a1) <- c("Chr","WinCenter","Pi","Intro","Win","Peaks")
colnames(a2) <- c("Chr","WinCenter","Pi","Intro","Win","Peaks")
colnames(n) <- c("Chr","WinCenter","Pi","Intro","Win","Peaks")
a1$Lineage <- "AME1"
a2$Lineage <- "AME2"
n$Lineage <- "NOV"
Df <- rbind(a1,a2,n)
Df$Class <- ifelse(Df$Peaks != "No-peak","High","Low")


p <- ggplot(Df) +
  geom_point(data = Df, aes(x = WinCenter, y = Pi, colour = Class),size=0.5) +
  scale_colour_manual(values = c("black","grey")) +
  scale_y_continuous(limits = c(0,0.018), breaks = c(0,0.01),labels = c("0.00","0.01")) +
  facet_grid(Lineage ~ Chr, scales = "free_x", shrink=FALSE) +
  scale_x_continuous(breaks = c(0,2000000,4000000,6000000),labels = c(0,2,4,6)) +
  labs(y = "Nucl. diversity") +
  theme(panel.background = element_rect(fill="grey95"),
        panel.grid = element_blank(),
        axis.text = element_text(size=10),
        axis.title.y = element_text(size=14),
        axis.title.x = element_blank(),
        strip.text.x = element_blank(),
        strip.text.y = element_text(angle = 0, size = 12),
        strip.background = element_rect(colour=NA, fill="white"),
        panel.border = element_blank(),axis.line.x = element_line(),axis.line.y = element_line(),
        panel.spacing.x = unit(1,"line"),
        panel.spacing.y = unit(0.1,"line"),
        legend.position = "None")
p
gt = ggplotGrob(p)
#gt$widths
gt$widths[5] = unit(1, "null")
gt$widths[7] = unit(0.98, "null")
gt$widths[9] = unit(0.53, "null")
gt$widths[11] = unit(0.49, "null")
gt$widths[13] = unit(0.41, "null")
gt$widths[15] = unit(0.4, "null")
gt$widths[17] = unit(0.4, "null")
gt$widths[19] = unit(0.36, "null")
grid.arrange(gt)


### Plotting ULM ancestry

anc <- read.table("createTable.txt",sep = "\t",header = TRUE)

a <- ggplot(anc) + aes(x = Start, y = Freq) + 
  scale_y_continuous(limits = c(0,0.5), breaks = c(0.5),labels = c("0.50")) +
  scale_x_continuous(breaks = c(0,2000000,4000000,6000000),labels = c(0,2,4,6)) +
  geom_area(fill = "#1874CD") +
  facet_grid(Lineage ~ Chrom, scales = "free_x", shrink = FALSE) +
  labs(y = "ULM ancestry") +
  theme(panel.background = element_rect(fill="grey95"),
        panel.grid = element_blank(),
        axis.text = element_text(size=10),
        axis.title.y = element_text(size=14),
        axis.title.x = element_blank(),
        strip.text.x = element_blank(),
        strip.text.y = element_text(angle = 0, size = 12),
        strip.background = element_rect(colour=NA, fill="white"),
        panel.border = element_blank(),axis.line.x = element_line(),axis.line.y = element_line(),
        panel.spacing.x = unit(1,"line"),
        panel.spacing.y = unit(0.1,"line"),
        legend.position = "None")
a
ga = ggplotGrob(a)
ga$widths
ga$widths[5] = unit(1, "null")
ga$widths[7] = unit(0.98, "null")
ga$widths[9] = unit(0.53, "null")
ga$widths[11] = unit(0.49, "null")
ga$widths[13] = unit(0.41, "null")
ga$widths[15] = unit(0.4, "null")
ga$widths[17] = unit(0.4, "null")
ga$widths[19] = unit(0.36, "null")
grid.arrange(ga)



### Getting and plotting MK targets

mk <- read.table("./createBed_targets_count.bed",sep="\t",header=FALSE)
colnames(mk) <- c('Chrom','win.start','win.stop','count')
mk$Start <- mk$win.start+1
mk$Lineage <- "AME2"  # name is just to adjust the space right of the plot

c <-ggplot(mk) + aes(x = Start, y = count) +
  geom_area(fill = "#D55E00") +
  facet_grid(Lineage ~ Chrom, scales = "free_x", shrink = FALSE) +
  scale_x_continuous(breaks = c(0,2000000,4000000,6000000),labels = c(0,2,4,6)) +
  scale_y_continuous(limits = c(0,8),breaks = c(0,4,8), labels = c("0.00","4.00","8.00")) +
  #scale_y_continuous(limits = c(0,8),breaks = c(0,4,8), labels = c("    0","    4","    8")) +
  labs(y = "MK") +
  theme(panel.background = element_rect(fill="grey95"),
        panel.grid = element_blank(),
        axis.text = element_text(size=10),
        axis.title.y = element_text(size=14),
        axis.title.x = element_blank(),
        strip.text.x = element_blank(),
        strip.text.y = element_text(angle = 0, size = 12,  color="white"),
        strip.background = element_rect(colour=NA, fill="white"),
        panel.border = element_blank(),axis.line.x = element_line(),axis.line.y = element_line(),
        #panel.border = element_rect(fill=NA),
        panel.spacing.x = unit(1,"line"),
        panel.spacing.y = unit(0.1,"line"),
        legend.position = "None")
c

gc = ggplotGrob(c)
#gc$widths
gc$widths[5] = unit(1, "null")
gc$widths[7] = unit(0.98, "null")
gc$widths[9] = unit(0.53, "null")
gc$widths[11] = unit(0.49, "null")
gc$widths[13] = unit(0.41, "null")
gc$widths[15] = unit(0.4, "null")
gc$widths[17] = unit(0.4, "null")
gc$widths[19] = unit(0.36, "null")
grid.arrange(gc)


### Getting and plotting Dxy

dxy <- read.csv("./combineTables.csv",header=TRUE)
dxy1 <- dxy %>% dplyr::select('gene_ID','Chrom','Gff_Start','dxy.ame1.ame2','high1.a1a2')
colnames(dxy1) <- c('gene_ID','Chrom','Gff_Start','Dxy','Status')
dxy1$Pair <- c("A1A2")
dxy2 <- dxy %>% dplyr::select('gene_ID','Chrom','Gff_Start','dxy.ame.nov','high1.anov')
colnames(dxy2) <- c('gene_ID','Chrom','Gff_Start','Dxy','Status')
dxy2$Pair <- c("ANOV")
dxy.m <- rbind(dxy1,dxy2)

x <- ggplot(dxy.m) +
  geom_point(data = dxy.m, aes(x = Gff_Start, y = Dxy, colour = Status),size=0.5) +
  scale_colour_manual(values = c("grey","black")) +
  scale_y_continuous(limits = c(0,0.04), breaks = c(0,0.03),labels = c("0.00","0.03")) +
  facet_grid(Pair ~ Chrom, scales = "free_x", shrink=FALSE) +
  scale_x_continuous(breaks = c(0,2000000,4000000,6000000),labels = c(0,2,4,6)) +
  labs(y = "Dxy") +
  theme(panel.background = element_rect(fill="grey95"),
        panel.grid = element_blank(),
        axis.text = element_text(size=10),
        axis.title.y = element_text(size=14),
        axis.title.x = element_blank(),
        strip.text.x = element_blank(),
        strip.text.y = element_text(angle=0,size = 12,colour="black"),
        strip.background = element_rect(colour=NA, fill="white"),
        panel.border = element_blank(),axis.line.x = element_line(),axis.line.y = element_line(),
        panel.spacing.x = unit(1,"line"),
        panel.spacing.y = unit(0.1,"line"),
        legend.position = "None")
x
gx = ggplotGrob(x)
#gt$widths
gx$widths[5] = unit(1, "null")
gx$widths[7] = unit(0.98, "null")
gx$widths[9] = unit(0.53, "null")
gx$widths[11] = unit(0.49, "null")
gx$widths[13] = unit(0.41, "null")
gx$widths[15] = unit(0.4, "null")
gx$widths[17] = unit(0.4, "null")
gx$widths[19] = unit(0.36, "null")
grid.arrange(gx)




### Plotting all graphs together


pdf('plot_TargetsAndDiversity.pdf',w=10,h=7)
lay <- rbind(c(1),c(1),c(1),c(1),c(1),c(2),c(2),c(2),c(2),c(2),c(3),c(3),c(3),c(3),c(3),c(4),c(4))
grid.arrange(ga,gt,gx,gc,ncol = 1, layout_matrix = lay,
             bottom = "Chromosome position (Mb)")
dev.off()

