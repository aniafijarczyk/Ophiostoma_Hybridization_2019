library(ggplot2)
library(gridExtra)
library(dplyr)



### Getting strain order
samp <- read.table(file = "strains_lineage_order.txt", sep = "\t", header = TRUE)
sample_order <- as.vector(samp$StrainOrder)

### Getting table with mating type tables
d <- read.table("MAT_loci_depth.txt",sep='\t',header=TRUE)
v <- d[c(6,8,10,12,14,15,17,18,19,20,21)]
m <- v %>% mutate(mat112=meanDepth.MAT112/Tot_mean_depth,mat113=meanDepth.MAT113/Tot_mean_depth,mat111=meanDepth.MAT111/Tot_mean_depth,
                  mat12=meanDepth.MAT12/Tot_mean_depth)
g <- m %>% gather(locus, depth, mat112:mat111)


### Plotting coverage depths
p1 <- ggplot(g) +
  aes(locus,name.pop) +
  geom_tile(aes(fill=depth),color="white") +
  scale_fill_gradient(low = "white", high = "darkblue") +
  scale_y_discrete(limits = rev(sample_order)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.background = element_blank())+
  scale_x_discrete(position = "top",labels=c("MAT1-1-1  ","MAT1-1-2  ","MAT1-1-3  "))+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_text(size=10),
        axis.text.x=element_text(size=11),
        legend.position = "bottom",
        legend.key.height=unit(0.3,"cm"),
        axis.ticks = element_blank())
p1

### Plotting mating type
p2 <- ggplot(g) + 
  aes(MThead,name.pop) +
  geom_tile(aes(fill=Mating.Type),color="white") +
  scale_fill_manual(values = c("MAT1" = "deeppink3", "MAT2" = "navy")) +
  scale_y_discrete(limits = rev(sample_order)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.background = element_blank())+
  scale_x_discrete(position = "top",labels = c("Mating type"))+
  #geom_hline(yintercept=c(21.5,48.5,75.5))+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_text(size=11),
        legend.position = "bottom",
        legend.title = element_blank(),
        axis.ticks = element_blank())
p2

### Plotting structure
st <- read.table(file = "collectWindow_changedNames.out", header = FALSE, sep = "\t", col.names = c("Sample","Chrom","Str_nov","Str_ulm"))
stg <- st %>% gather(anc = Str_nov:Str_ulm, ancestry)
stg$Sample <- factor(stg$Sample, levels = rev(sample_order))

p3 <- ggplot(stg) + 
  aes(x = Sample, y = value, fill = ancestry) +
  geom_bar(position = "fill",stat = "identity") +
  scale_fill_manual(values = c("khaki3", "#1874CD"), labels = c("NOV","ULM")) +
  scale_x_discrete(position = "top") +
  scale_y_continuous(position = "bottom", labels = c("","","Admixture","","")) +
  coord_flip() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_text(size=11),
        legend.position = "bottom",
        legend.title = element_blank(),
        axis.ticks = element_blank())
p3

pdf(file="plotDepth_heatmap.pdf",width = 6.5, height = 16)  
grid.arrange(p1,p3,p2,nrow=1,widths=c(0.3,0.12,0.12),heights=c(1))
dev.off()

