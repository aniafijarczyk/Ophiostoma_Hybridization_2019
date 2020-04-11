library(ggplot2)
library(dplyr)
library(viridis)
library(gridExtra)




# Setting colours

pal16 <- c("#9986A5","#CB2314","#899DA4","#00A08A","#1E1E1E","#DC863B","#EABE94",
           "#0B775E","#9A8822","#F8AFA8","#FD6467","#C6CDF7","#85D4E3","#046C9A",
           "#FDD262","#FAD510")



# Reading files

pops <- read.table("./data/Species_distinction_PCA.txt", sep = "\t", header = TRUE, na.strings = c("NaN"))
popGeo <- pops[c(14,25)]
df <- read.table("./data/selectStrainsForbranchEstimationBootstrapAllOU.out", sep = "\t", header = TRUE)

regs <- read.table("./data/regions_dataframe.txt",sep = "\t",header = TRUE)
regs["mid"] <- (regs["stop"]-regs["start"]+1)/2 + regs["start"]



# Preparing data frame

strains <- df$strain_x
locs <- df$region
sl <- c()
for (i in 1:length(strains)) {
  x = paste0(strains[i],"_",locs[i])
  sl[i] <- x
}
df$strain_loc <- sl
df["Length"] <- df["tot"]

dg <- df %>% group_by(strain_loc,strain_x,Length) %>% summarize(mean = mean(ttout))
dm <- merge(df,dg,by = c("strain_loc","Length"), sort = FALSE)
dm["diff"] <- dm["mean"] - dm["ttout"]
db <- dm %>% group_by(strain_loc,mean,Length,chrom,region,strain_x.x)%>% summarize(p5=quantile(diff,probs=0.05),p95=quantile(diff,probs=0.95))
db["low"] <- db["mean"] - db["p95"]
db["high"] <- db["mean"] - db["p5"]



# Plotting bootstrap distribution by region

dsort <- dm[order(dm$mean),] 
sordered <- dsort$strain_loc
ordered <- unique(sordered)
dm$strain_loc <- factor(dm$strain_loc, levels = ordered)

p1 <- ggplot(dm) +
  geom_jitter(data = dm, aes(x = strain_loc, y = ttout, colour = as.factor(region)), width = 0.3, size = 0.5) +
  scale_color_manual(values = pal16) +
  labs(x = "Strain_location", y = "relative branch length", colour = "Region") +
  geom_errorbar(data = db, aes(x = strain_loc, ymin = low, ymax = high), width = 0.3) + 
  geom_point(data = db, aes(x = strain_loc, y = mean), size = 1.5, shape = 23, fill = "black", colour = "black") +
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        panel.background = element_blank(),
        axis.line = element_line(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 12),
        legend.key=element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=5)))
p1

png("plot_IR_length_by_region.png",w=9,h=5,res=150,units="in")
p1
dev.off()



# Plot regions on chromsoosmes

chrs <- data.frame("chr"=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8"),
                   "start"=c(1,1,1,1,1,1,1,1),
                   "stop"=c(6937932,6817711,3669772,3419703,2848703,2801594,2758224,2531247),
                   "stat"=c(0,0,0,0,0,0,0,0))

p2 <- ggplot(chrs) + 
  geom_rect(data = chrs, aes(xmin=start,xmax=stop,ymin=-1,ymax=1), colour="black", fill="white")+
  facet_wrap(~chr,ncol = 1, strip.position = "left")+
  geom_rect(data = regs, aes(xmin=start,xmax=stop,ymin=-1,ymax=1,group = chr, fill = as.factor(region))) +
  scale_fill_manual(values = pal16) +
  geom_text(data = regs, aes(x = mid, y = 0, group = chr, label = as.factor(region)), colour = c(rep("black",11),"white",rep("black",2)))+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = "None",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        strip.background = element_blank(),
        strip.text.y = element_text(angle = 180, size = 14))
p2

png("plot_regions.png",w=6,h=4,res=150,units="in")
p2
dev.off()

