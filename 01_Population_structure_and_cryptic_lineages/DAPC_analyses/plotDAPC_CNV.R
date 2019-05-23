
library(vcfR)
library(adegenet)


### Getting vcf
vcf <- read.vcfR("filtered_0.2-himal.gt.vcf", verbose = FALSE)
x1 <- vcfR2genlight(vcf)

### Running dapc
grp1=find.clusters(x1,max.n.clust = 40, parallel=FALSE)
60
5

dapc1=dapc(x1,grp1$grp, parallel=FALSE)
60
5

head(dapc1)
coord <- dapc1$ind %>% data.frame()
coord$ID <- row.names(dapc1$ind)
#write.table(coord, "coord_DAPC")


### Getting population info
IdInfo <- read.table("Species_distinction_PCA.txt", header = T) %>% mutate(ID = sub(".*OHT_[0-9]+?.(.+?)_aln.*", "\\1", SpeciesID), StrainGrp = paste(ID, Group, sep = "_"), Idnb = Ind) %>% select(-SpeciesID, -PCA1, -PCA2)
head(IdInfo)
coord_info <- left_join(coord, IdInfo, by = "ID")
head(coord_info)

dapc1$eig
dapc1$eig[1]/sum(dapc1$eig)
dapc1$eig[2]/sum(dapc1$eig)
dapc1$eig[3]/sum(dapc1$eig)

dapc1$pca.eig
dapc1$pca.eig[1]/sum(dapc1$pca.eig)
dapc1$pca.eig[2]/sum(dapc1$pca.eig)
dapc1$pca.eig[3]/sum(dapc1$pca.eig)



### Getting lineage names

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
sort <- coord_info[order(coord_info[,30]),]


### Plotting DAPC

PC_1_2 <- ggplot(sort, aes( x = LD1.x, y = LD2.x, col = groupnames)) + 
  geom_vline(xintercept = 0, cex = 0.2) + 
  geom_hline(yintercept = 0, cex = 0.2) + 
  geom_point(alpha = 0.8, cex = 5) + 
  theme_bw() + 
  scale_color_manual(values = c(cl[136], cl[149], cl[497], cl[131])) + 
  xlab("DPC1 (92.5%)") + 
  ylab("DPC2 (3.6%)") +
  theme(panel.grid = element_blank(),
        legend.title = element_blank())
PC_1_2


PC_2_3 <- ggplot(sort, aes( x = LD2.x, y = LD3, col = groupnames)) + 
  geom_vline(xintercept = 0, cex = 0.2) + 
  geom_hline(yintercept = 0, cex = 0.2) + 
  geom_point(alpha = 0.8, cex = 5)+ 
  theme_bw() + 
  scale_color_manual(values = c(cl[136], cl[149], cl[497], cl[131])) + 
  xlab("DPC2 (3.6%)") + 
  ylab("DPC3 (2.4%)") +
  theme(panel.grid = element_blank(),
        legend.title = element_blank())
PC_2_3

pdf("DAPC_CNV_DPCplots.pdf", width = 4.5, height = 6)
grid.arrange(PC_1_2, PC_2_3)
dev.off()


### Plotting eigenvalues and cumulative distribution of DPCs

PC_eig <- cbind("eig" = dapc1$pca.eig, "PC" = c(1:length(dapc1$pca.eig))) %>% as.data.frame() %>% mutate(PC_statut = ifelse (PC > 60, "retain", "discard"), var = eig/sum(eig)) %>% mutate(cum_sum = cumsum(var))

PC_eig_plot <- ggplot(PC_eig, aes(x = PC, y = cum_sum, fill = PC_statut)) + 
  geom_bar(stat = "identity") + 
  theme_bw() + 
  scale_fill_manual(values = c("grey30", "grey")) + 
  theme(legend.position = "none") + 
  ylab("Cumulative sum of variance") +
  theme(panel.grid = element_blank())
PC_eig_plot

DPC_eig <- cbind("eig" = dapc1$eig, "DPC" = c(1:length(dapc1$eig))) %>% as.data.frame() %>% mutate(var = eig/sum(eig)*100) %>% mutate(cum_sum = cumsum(var))
DPC_eig_plot <- ggplot(DPC_eig, aes(x = DPC, y = var)) + 
  geom_bar(stat = "identity") + 
  theme_bw() + 
  theme(legend.position = "none") + 
  ylab("Percentage of variance") +
  theme(panel.grid = element_blank())
DPC_eig_plot



pdf("DAPC_CNV_eigenvalues.pdf", width = 5, height = 6)
grid.arrange(DPC_eig_plot, PC_eig_plot)
dev.off()





