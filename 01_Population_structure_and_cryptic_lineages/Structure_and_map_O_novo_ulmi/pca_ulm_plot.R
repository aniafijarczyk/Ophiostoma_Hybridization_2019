library(vcfR)
library(ggplot2)
library(adegenet)
library(ggrepel)


vcf <- read.vcfR("SNP_5000_ulmi.vcf", verbose = FALSE)
x <- vcfR2genlight(vcf)

infos=read.table('strains_pops.txt',sep='\t',header=TRUE)

x@ploidy=rep(as.integer(1),21)
pca=glPca(x, parallel = F)
2
scatter(pca)
scores=as.data.frame(pca$scores)
scores = tibble::rownames_to_column(scores, "ID")


d=merge(scores,infos,by.x="ID",by.y="NameBam",all.x=T,all.y=F,sort=FALSE)

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
