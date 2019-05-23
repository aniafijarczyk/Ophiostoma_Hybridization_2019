library(dplyr)
library(ggplot2)
library(gridExtra)


### Importing structure results

k2 <- read.csv("./Structure_genomewide_results_novoulmi_K2.csv",sep=",",header=TRUE)
k3 <- read.csv("./Structure_genomewide_results_novoulmi_K3.csv",sep=",",header=TRUE)
k4 <- read.csv("./Structure_genomewide_results_novoulmi_K4.csv",sep=",",header=TRUE)
k5 <- read.csv("./Structure_genomewide_results_novoulmi_K5.csv",sep=",",header=TRUE)
pop <- read.table("strains_pops.txt",sep="\t",header=TRUE)

k2pop <- merge(k2,pop,by.x = 'Label',by.y = 'SpName')
k3pop <- merge(k3,pop,by.x = 'Label', by.y = 'SpName')
k4pop <- merge(k4,pop,by.x = 'Label', by.y = 'SpName')
k5pop <- merge(k5,pop,by.x = 'Label', by.y = 'SpName')


### Encoding regions

k5pop$Region.x <- as.character(k5pop$Region.x)
k5pop[28,10] <- 'Unknown'
k5pop$Region.x <- as.factor(k5pop$Region.x)

nov <- k5pop %>% filter(Blank == 0)
regions <- nov$Region.x
regcode <- vector("integer", length(regions))
for (i in 1:length(regions)) {
  if (regions[i] == 'Namerica') {
    x <- 1
  } else if (regions[i] == 'Europe') {
    x <- 2
  } else if (regions[i] == 'Asia') {
    x <- 3
  } else if (regions[i] == 'NewZealand') {
    x <- 4
  } else {x <- 5}
  regcode[i] <- x
}
nov$regcode <- regcode

### Sorting data frame

sortedDF <- nov[order(nov[,23], nov[,6], nov[,4], nov[,5], nov[,7], nov[,8]),]
o <- sortedDF['NameBam']
colnames(o) <- 'row_order'

### Setting colors

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
f2 <- dk2 %>% filter(Blank == 0)
dn <- f2[c('row_order','AME','NOV')]
gk2 <- dn %>% gather(key = 'Species', value = 'Admixture', c(2,3))
# reorder dataframe
gk2$row_order <- factor(gk2$row_order, levels = as.vector(o$row_order))

p2 <- ggplot(gk2) + aes(x = row_order, y = Admixture, fill = Species) +
  geom_bar(stat = 'identity', col=NA, width=1) +
  scale_fill_manual(values = c(olivedrab4,orange3)) +
  labs(y = "K = 2") +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none",
        #axis.text.x = element_text(angle=90),
        axis.text.x = element_blank())
p2



### K3 plot

dk3 <- merge(o, k3pop, by.x = 'row_order', by.y = 'NameBam')
f3 <- dk3 %>% filter(Blank == 0)
dn <- f3[c('row_order','AME1','AME2','NOV')]
gk3 <- dn %>% gather(key = 'Species', value = 'Admixture', c(2,3,4))
# reorder dataframe
gk3$row_order <- factor(gk3$row_order, levels = as.vector(o$row_order))

p3 <- ggplot(gk3) + aes(x = row_order, y = Admixture, fill = Species) +
  geom_bar(stat = 'identity', col=NA, width=1) +
  scale_fill_manual(values = c(firebrick3,goldenrod2,olivedrab4)) +
  labs(y = "K = 3") +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none",
        #axis.text.x = element_text(angle=90),
        axis.text.x = element_blank())
p3


### K4 plot

dk4 <- merge(o, k4pop, by.x = 'row_order', by.y = 'NameBam')
f4 <- dk4 %>% filter(Blank == 0)
dn <- f4[c('row_order','C1','C2','C3','C4')]
gk4 <- dn %>% gather(key = 'Species', value = 'Admixture', c(2,3,4,5))
# reorder dataframe
gk4$row_order <- factor(gk4$row_order, levels = as.vector(o$row_order))

p4 <- ggplot(gk4) + aes(x = row_order, y = Admixture, fill = Species) +
  geom_bar(stat = 'identity', col=NA, width=1) +
  scale_fill_manual(values = c('aquamarine3',firebrick3,olivedrab4,goldenrod2)) +
  labs(y = "K = 4") +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none",
        #axis.text.x = element_text(angle=90),
        axis.text.x = element_blank())
p4



### K5 plot

dk5 <- merge(o, k5pop, by.x = 'row_order', by.y = 'NameBam')
f5 <- dk5 %>% filter(Blank == 0)
dn <- f5[c('row_order','C1','C2','C3','C4','C5')]
gk5 <- dn %>% gather(key = 'Species', value = 'Admixture', c(2,3,4,5,6))
# reorder dataframe
gk5$row_order <- factor(gk5$row_order, levels = as.vector(o$row_order))

p5 <- ggplot(gk5) + aes(x = row_order, y = Admixture, fill = Species) +
  geom_bar(stat = 'identity', col=NA, width=1) +
  scale_fill_manual(values = c('aquamarine3','darkcyan',olivedrab4,goldenrod2,firebrick3)) +
  labs(y = "K = 5") +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none",
        axis.text.x = element_blank())
p5



### Plotting regions

st <- sortedDF %>% filter(Blank == 0)
st$Stat <- 1.0
st$SamName <- factor(st$Name, levels = as.vector(o$row_order))

pr <- ggplot(st) + aes(x = SamName, y = Stat, fill = as.factor(regcode)) +
  geom_bar(stat = 'identity', col=NA, width=1) +
  geom_text(x = 18.5, y = 0.5, label = "North America") + 
  geom_text(x = 51.5, y = 0.5, label = "Europe") + 
  geom_text(x = 66.5, y = 0.5, label = "Asia",colour = "white") + 
  geom_text(x = 71.5, y = 0.5, label = "NZ",colour = "white",angle=90) +
  geom_text(x = 73, y = 0.5, label = "unk",angle=90) +
  scale_fill_manual(values = c("grey80","grey60","grey40","grey20","grey90")) +
  labs(y = "R") +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(colour = "white"),
        legend.position = "none",
        #axis.text.x = element_text(angle=90),
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour = "white"))
pr

### Plotting lineages

pl <- ggplot(st) + aes(x = SamName, y = Stat, fill = GroupSamtools) +
  geom_bar(stat = 'identity', col=NA, width=1) +
  scale_fill_manual(values = c(firebrick3,goldenrod2,olivedrab4)) +
  scale_x_discrete(labels = st$Label) +
  labs(y = "L") +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(colour = "white"),
        legend.position = "none",
        axis.text.x = element_text(angle=90),
        axis.text.y = element_text(colour = "white"))
pl




### Plotting all graphs together

grid.arrange(p2,p3,p4,p5,pr,pl,ncol=1)

pdf(file="Figure_structure_regions.pdf",width=12,height=8)
lay <- rbind(c(1),c(1),c(2),c(2),c(3),c(3),c(4),c(4),c(5),c(6),c(6))
grid.arrange(p2,p3,p4,p5,pr,pl, ncol = 1, layout_matrix = lay)
dev.off()
