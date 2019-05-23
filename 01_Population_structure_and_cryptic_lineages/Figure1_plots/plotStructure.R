library(dplyr)
library(ggplot2)
library(tidyr)
library(reshape)
library(gridExtra)

### Import structure results

k2 <- read.csv("Structure_genomewide_results_K2.csv",sep=",",header=TRUE)
k3 <- read.csv("Structure_genomewide_results_novoulmi_K3.csv",sep=",",header=TRUE)
k2n <- read.csv("Structure_genomewide_results_novoulmi_K2.csv",sep=",",header=TRUE)
pop <- read.table("Species_distinction_PCA.txt", sep="\t", header = T)

k2pop <- merge(k2,pop,by.x = 'Label',by.y = 'SpName')
k3pop <- merge(k3,pop,by.x = 'Label', by.y = 'SpName')
k2npop <- merge(k2n,pop,by.x = 'Label', by.y = 'SpName')
  
### Import  order
o <- read.table("plotTree_strain_order_corr.txt",sep="\t",header=FALSE)
colnames(o) <- 'row_order'


### Colors
palegreen <- colors()[517] 
goldenrod2 <- colors()[149]
olivedrab4 <- colors()[497]
firebrick3 <- colors()[136]
aquamarine4 <- colors()[12]
dodgerblue3 <- colors()[131]
orange3 <- colors()[501]
khaki3 <- colors()[385]

### K2 plot

dk2 <- merge(o, k2pop, by.x = 'row_order', by.y = 'NameBam')
dn <- dk2[c('row_order','ULM','Novo.ulmi')]
gk2 <- dn %>% gather(key = 'Species', value = 'Admixture', c(2,3))
# reorder dataframe
gk2$row_order <- factor(gk2$row_order, levels = as.vector(o$row_order))

p1 <- ggplot(gk2) + aes(x = row_order, y = Admixture, fill = Species) +
  geom_bar(stat = 'identity', col=NA, width=1) +
  scale_fill_manual(values = c(khaki3,dodgerblue3)) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.position = "none",
        #axis.text.x = element_text(angle=90),
        axis.text = element_blank())
p1

### K2 nov plot

dk2n <- merge(o, k2npop, by.x = 'row_order', by.y = 'NameBam')
dnn <- dk2n[c('row_order','AME','NOV','Blank')]
gk2n <- dnn %>% gather(key = 'Species', value = 'Admixture', c(2,3))
# reorder dataframe
gk2n$row_order <- factor(gk2n$row_order, levels = as.vector(o$row_order))

p3 <- ggplot(gk2n) + aes(x = row_order, y = Admixture, fill = Species) +
  geom_bar(stat = 'identity', col=NA, width=1) +
  scale_fill_manual(values = c(olivedrab4,orange3)) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.position = "none",
        #axis.text.x = element_text(angle=90),
        axis.text = element_blank())
p3



### K3 plot

dk3 <- merge(o, k3pop, by.x = 'row_order', by.y = 'NameBam')
dc <- dk3[c('row_order','AME1','AME2','NOV','Blank')]
gk3 <- dc %>% gather(key = 'Species', value = 'Admixture', c(2,3,4,5))
gk3$row_order <- factor(gk3$row_order, levels = as.vector(o$row_order))

p2 <- ggplot(gk3) + aes(x = row_order, y = Admixture, fill = Species) +
  geom_bar(stat = 'identity', col=NA, width=1) +
  scale_fill_manual(values = c("#CD2626","#EEB422","white","#698B22")) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.position = "none",
        axis.text = element_blank())
p2




### Combine the plots

pdf("plotStructure.pdf", w = 8, h = 4)
grid.arrange(p2,p3,p1,ncol=1)
dev.off()


