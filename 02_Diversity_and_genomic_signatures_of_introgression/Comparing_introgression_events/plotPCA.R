rm(list=ls())
library("SNPRelate")
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(gridExtra)
library("wesanderson")



# Setting colour palette

pal1 <- wes_palette("Darjeeling1")



# Reading vcf files (example for location 15)

loc=15
chrloc = as.character(loc)
vcf.fn <- paste0("./data/location_",chrloc,".recode.vcf")
snpgdsVCF2GDS_R(vcf.fn, paste0("./data/location_",chrloc,".recode.gds"), method="biallelic.only",option=snpgdsOption(OphioH327chr_1=1, OphioH327chr_2=2, OphioH327chr_3=3, OphioH327chr_4=4, OphioH327chr_5=5, OphioH327chr_6=6, OphioH327chr_7=7, OphioH327chr_8=8))
genofile <- snpgdsOpen(paste0("./data/location_",chrloc,".recode.gds"))



# Running PCA

pca <- snpgdsPCA(genofile, num.thread=3)
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))
tab <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)
head(tab)
plot(tab$EV2, tab$EV1, xlab="eigenvector 2", ylab="eigenvector 1")



# Getting sample ids and population info

sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
dm <- read.table("./data/Species_distinction_PCA.txt",sep="\t",header=TRUE)
pop_code <- dm[c(3,8,12,20,21)]
colnames(pop_code) <- c("sample.id","Strain","GeoGroup","Lineage","SamName")
pop_code$Geo <- ifelse(pop_code$GeoGroup == "A", "NAmerica","Europe")
pop_code$Temp <- ifelse(pop_code$Lineage == "U","ulmi","novo-ulmi")
geogroup <- pop_code$Geo
lines <- pop_code$Lineage
gr <- c()
for (i in 1:length(geogroup)) {
  x = paste0(lines[i],"_",geogroup[i])
  gr[i] <- x
}
pop_code$Group <- gr
pop_codev <- as.vector(pop_code$pop_code)
country <- as.vector(pop_code$GeoGroup)
head(cbind(sample.id, pop_codev))

tab <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  EV3 = pca$eigenvect[,3],
                  EV4 = pca$eigenvect[,4],
                  EV5 = pca$eigenvect[,5],
                  EV6 = pca$eigenvect[,6],
                  stringsAsFactors = FALSE)
head(tab)
tabp <- merge(tab, pop_code, by = 'sample.id', sort = FALSE)
head(tabp)



# Plotting PCA

cols1 <- c("#E69F00","#009E73","#56B4E9","#0072B2","#D55E00","#CC79A7","#E69F00")
cols2 <- c("#E69F00","#D55E00","#009E73","#56B4E9","#0072B2","#D55E00","#CC79A7","#E69F00")
pc1 <- as.character(round(pc.percent[1],2))
pc2 <- as.character(round(pc.percent[2],2))

# remember to adjust title of the plot
titles <- c("chr1:0-200000[5]","chr2:5400000-5600000","chr2:6200000-6600000[2]","chr7:3600000-3800000[2]","chr6:400000-600000[4]","chr8:600000-800000[5]")

# ... for all locations except 4

#ggplot(tabp) + aes(x = EV1, y = EV2, color = Group, label = Strain,shape=Group,alpha=Group) +
p <- ggplot(tabp) + aes(x = EV1, y = EV2, color = Group, shape=Group,alpha=Group) +
  geom_point(size=7,stroke=1) +
  scale_shape_manual(values=c(3,19,19)) +
  scale_alpha_manual(values=c(1,0.75,0.75)) +
  scale_colour_manual(values = c("black","#56B4E9","#0072B2")) +
  #geom_text(hjust=0, vjust=0) +
  labs(x =paste0("PC1=",pc1,"%"), y = paste0("PC1=",pc2,"%")) +
  ggtitle(titles[6]) + 
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(fill = FALSE),
        axis.title = element_text(size=14),
        axis.text = element_text(size=14),
        plot.title = element_text(size=18))
p

png(paste0("plot_loc",chrloc,".png"),w=9,h=5,res=150,units="in")
p
dev.off()



# Plotting PCA for location 4

#ggplot(tabp) + aes(x = EV1, y = EV2, color = Group, label=Strain,shape=Group,alpha=Group) +
ggplot(tabp) + aes(x = EV1, y = EV2, color = Group, shape=Group,alpha=Group) +
  geom_point(size=7,stroke=1) +
  scale_shape_manual(values=c(4,2,3,19,19)) +
  scale_alpha_manual(values=c(1,1,1,0.75,0.75)) +
  scale_colour_manual(values = c("#E69F00","#D55E00","black","#56B4E9","#0072B2")) +
  #geom_text(hjust=0, vjust=0) +
  labs(x =paste0("PC1=",pc1,"%"), y = paste0("PC1=",pc2,"%")) +
  ggtitle(titles[2]) + 
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(fill = FALSE),
        axis.title = element_text(size=14),
        axis.text = element_text(size=14),
        plot.title = element_text(size=18))

