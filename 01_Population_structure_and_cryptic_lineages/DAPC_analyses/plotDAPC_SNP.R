library(dplyr)
library(ggplot2)


### Reading files with DAPC coordinates
coord <- read.table("Ind_coord_dapcsamtools4.txt",sep=" ",header=F,skip=1) # DAPC for O. novo-ulmi only
colnames(coord) <- c("ID","DPC1","DPC2")
coord$ID <- as.factor(coord$ID)


### Reading population info
IdInfo <- read.table("Species_distinction_PCA.txt", header = T) %>% mutate(ID = sub(".*OHT_[0-9]+?.(.+?)_aln.*", "\\1", SpeciesID), StrainGrp = paste(ID, Group, sep = "_"), Idnb = Ind) %>% select(-SpeciesID, -PCA1, -PCA2)
IdInfo$SpeciesID.1 <- as.factor(IdInfo$SpeciesID.1)
head(IdInfo)

coord_info <- merge(coord, IdInfo, by.x = "ID", by.y = "NameBam")
head(coord_info)


### Getting population names
cl <- colors()[]
groups <- coord_info$GroupSamtools
grcode <- vector("integer", length(groups))
groupnames <- vector("integer", length(groups))
for (i in 1:length(groups)) {
  if (groups[i] == 'A2') {
    x <- 1
    y <- 'AME2'
  } else if (groups[i] == 'A1') {
    x <- 2
    y <- 'AME1'
  } else if (groups[i] == 'N') {
    x <- 3
    y <- 'NOV'
  } else if (groups[i] == 'U') {
    x <- 4
    y <- 'ULM'
  }
  grcode[i] <- x
  groupnames[i] <- y
}
coord_info$grcode <- grcode
coord_info$groupnames <- groupnames


### Plotting PCs

pc12 <- ggplot(coord_info, aes( x = DPC1, y = DPC2, col = groupnames)) + 
  geom_vline(xintercept = 0, cex = 0.2) + 
  geom_hline(yintercept = 0, cex = 0.2) + 
  geom_point(alpha = 0.8, cex = 5) + 
  theme_bw() + 
  scale_color_manual(values = c(cl[136], cl[149], cl[497], cl[131])) + 
  xlab("DPC1") + 
  ylab("DPC2") +
  theme(panel.grid = element_blank(),
        legend.title = element_blank())
pc12


pdf("DAPC_SNP_plots.pdf", width = 4, height = 4)
pc12
dev.off()

  