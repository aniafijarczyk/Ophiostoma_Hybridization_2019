#!/usr/bin/env RScript

# USAGE:
# ./runPCA.R
### ./runPCA.R raw_file=raw_file label_file=label_file pop_file=pop_file plot_order=plot_order output_dir=output_dir
#
# raw_file = "OphioH327chr_1_0a.raw"
# label_file = "Raw_label"
# pop_file = "Species_distinction_PCA.txt"
# plot_order = "plotLabelOrder_2.txt"
# output_dir="./"


argz <- commandArgs()
argz <- c("raw_file=OphioH327chr_1_0a.raw","label_file=Raw_label", "pop_file=Species_distinction_PCA.txt","plot_order=plotLabelOrder_2.txt","output_dir=./")

for(i in argz){
    if(strsplit(i, "=")[[1]][1]=="raw_file"){raw_file = strsplit(i, "=")[[1]][2]}
    if(strsplit(i, "=")[[1]][1]=="label_file"){label_file = strsplit(i, "=")[[1]][2]}
	if(strsplit(i, "=")[[1]][1]=="pop_file"){pop_file = strsplit(i, "=")[[1]][2]}
    if(strsplit(i, "=")[[1]][1]=="plot_order"){plot_order = strsplit(i, "=")[[1]][2]}
    if(strsplit(i, "=")[[1]][1]=="output_dir"){output_dir = strsplit(i, "=")[[1]][2]}
}



# Used libraries
library(adegenet)
library(dplyr)
library(ggplot2)



#getwd()

print(raw_file)
PlotLabel_order <- read.table(plot_order, header = F, col.names = "strain")
label <- read.table(label_file, header = T, col.names = "SpeciesID")
label$SpeciesID <- as.character(label$SpeciesID)
pop <- read.table(pop_file, header = T)
head(pop)

label_pop <- left_join(label, pop) %>% select(SpeciesID_samtools, Group.Samtools)
head(label_pop)

windows <- sub(".*/(Ophio.*).raw", "\\1", raw_file)

dt <- read.PLINK(raw_file, map.file = NULL, quiet = FALSE, chunkSize = 1000, parallel = require("parallel"), n.cores = 3, ploidy = 1)


#glPlot(dt)
pca0 <- glPca(dt)
2
#scatter(pca0)
pca0_df <- as.table(pca0$scores) %>% as.data.frame %>% reshape(idvar = "Var1", timevar = "Var2", direction = "wide")
colnames(pca0_df) <- c("SpeciesID_samtools", "PC1", "PC2")
pca0_pop <- left_join(pca0_df, pop) %>% select("SpeciesID_samtools", "PC1", "PC2", "Group.Samtools")

pdf(paste0(output_dir, windows, ".pdf"))
ggplot(pca0_pop, aes(x = PC1, y = PC2, col = Group.Samtools)) + geom_point(cex = 5, alpha = 0.5) + theme_bw() + scale_color_manual(values = c("firebrick3", "goldenrod2", "olivedrab4", "dodgerblue3")) + ggtitle(windows)
dev.off()

PC1_mean_ULM <- (pca0_pop %>% group_by(Group.Samtools) %>% summarise(PC1_mean = mean(PC1)) %>% as.data.frame %>% filter(Group.Samtools == "ULM"))[, 2]

pca0_dist <- pca0_pop %>% mutate(PC1_mean_ULM = PC1_mean_ULM, dist = abs(PC1_mean_ULM - PC1), file = sub(".*/(Ophio.*).raw", "\\1", raw_file), strain = sub(".*_[0-9]+?[.](.*)_calmd.*", "\\1", SpeciesID_samtools), dist_norm = dist / max(dist))

write.table(pca0_dist, file = paste0(output_dir, windows, ".txt"), quote = F, row.names = F)

if (file.exists(paste0(output_dir, "AllWindows.txt"))){
	write.table(pca0_dist, file = paste0(output_dir, "AllWindows.txt"), quote = F, row.names = F, col.names = F, append = T)
} else {
	write.table("SpeciesID PC1 PC2 Group.Samtools PC1_mean_ULM dist file strain dist_norm", file = paste0(output_dir,"AllWindows.txt"), quote = F, row.names = F, col.names = F)
	write.table(pca0_dist, file = paste0(output_dir, "AllWindows.txt"), quote = F, row.names = F, col.names = F, append = T)
}

#ggplot(pca0_dist, aes(x = file, y = strain)) + geom_raster(aes(fill = dist_norm)) + scale_fill_gradient(low = "dodgerblue3", high = "khaki3") + scale_y_discrete(limits = rev(PlotLabel_order$strain), expand = c(0,0))

