library(vcfR)
library(ape)
library(adegenet)
library(ggplot2)
library(dplyr)
library(readxl)
library(ggrepel)



vcf <- read.vcfR("SNP_5000_ulmi.vcf", verbose = FALSE)
x <- vcfR2genlight(vcf)

x@ploidy=rep(as.integer(1),21)

pca=glPca(x, parallel = F)
2
scatter(pca)

scores=as.data.frame(pca$scores)
scores = tibble::rownames_to_column(scores, "ID")

IDs=read_xlsx("IDs.xlsx")
infos= read.csv2("isolates_infos.csv")
infos=infos[,4:10]

d=merge(scores, IDs, by.x="ID", by.y="seq",all.x = T, all.y = F)
d=merge(d, infos, by.x="ID_simple",by.y="ID", all.x=T, all.y=F)

p <- ggplot(d, aes(as.numeric(as.character(PC1)), as.numeric(as.character(PC2))), label=country) + 
  geom_vline(xintercept=0, color = "grey80") + 
  geom_hline(yintercept=0, color = "grey80") + 
  geom_label_repel(aes(label = country),box.padding = 0.4, point.padding = 0.5,segment.color = 'grey50') + 
  geom_point(size = 4) + 
  theme_bw() + 
  ylab("PC2") + 
  xlab("PC1") +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 12),
        panel.grid = element_blank())

pdf(file="Figure_PCA_SNP_OU.pdf",width=5,height=5)
p
dev.off()
