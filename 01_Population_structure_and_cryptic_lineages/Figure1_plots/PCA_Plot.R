

library(ggplot2)
library(dplyr)
library(reshape)
library(gridExtra)



cl <- colors()[]
sample_order <- c("TR88", "R69", "H332", "H413a", "Yu141", "H346", "R108", "H323", "Yu16", "H371", "P122", "R22", "V19", "H986", "H260", "CKT4", "H294", "TR67", "H276", "H648", "P114", "V14", "H328", "BP14", "AT146", "AT83", "782401", "H340", "H312", "DDS81", "DDS93", "H2091", "H2095", "H618", "H623", "H234D", "NZFS601", "NZFS184C", "BP18", "DDS100", "ES130", "H295", "DDS82", "DDS154", "H2083", "DDS25", "DDS170", "DDS202", "H2086", "H239", "TB2", "TB5", "H2098", "H248")

# Reading files

scores_pca <- read.table("scores_pca.txt", header = T) %>% mutate(ID = sub(".*OHT_[0-9]+?.(.+?)_calmd.*", "\\1", SpeciesID.1))
#head(scores_pca)

IdInfo <- read.table("Species_distinction_PCA.txt", header = T) %>% mutate(ID = sub(".*OHT_[0-9]+?.(.+?)_aln.*", "\\1", SpeciesID), StrainGrp = paste(ID, Group, sep = "_"), Idnb = Ind) %>% select(-SpeciesID, -PCA1, -PCA2)
#head(IdInfo)

scores_pca_info <- left_join(scores_pca, IdInfo, by = "ID")
head(scores_pca_info)


# Plotting PCA

PC_1_2 <- ggplot(scores_pca_info, aes( x = PC1, y = PC2, col = GroupSamtools)) + geom_vline(xintercept = 0, cex = 0.2) + geom_hline(yintercept = 0, cex = 0.2) + geom_point(alpha = 0.5, cex = 5)+ theme_bw() + scale_color_manual(values = c(cl[501], cl[501], cl[497], cl[131]))

PC_2_3 <- ggplot(scores_pca_info, aes( x = PC2, y = PC3, col = GroupSamtools)) + geom_vline(xintercept = 0, cex = 0.2) + geom_hline(yintercept = 0, cex = 0.2) + geom_point(alpha = 0.5, cex = 5)+ theme_bw() + scale_color_manual(values = c(cl[136], cl[149], cl[497], cl[131]))

pdf("PCA_plots_ulmi_novsp.pdf", width = 5, height = 6)
grid.arrange(PC_1_2, PC_2_3)
dev.off()

