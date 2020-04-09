library(dplyr)
library(ggplot2)
library(tidyr)
library(reshape)
library(gridExtra)

### import structure results

k2 <- read.table("ulmi_k2.txt",sep="",skip = 1, header=FALSE)
k3 <- read.table("ulmi_k3.txt",sep="",skip = 1, header=FALSE)

pop <- read.table("strains_pops.txt",sep="\t",header=TRUE, na.strings = 'na')

k2pop <- merge(k2,pop,by.x = 'V2', by.y = 'NameBam')
k3pop <- merge(k3,pop,by.x = 'V2', by.y = 'NameBam')


### order strains

k2pop$Region <- as.factor(k2pop$Region)

df <- k2pop[order(k2pop[,17], k2pop[,5]),]
rownames(df) <- 1:nrow(df)
names <- df$SpName
continents <- df$Region
new_names <- paste(names,continents, sep="_") 
df['new_names'] <- new_names


### colors

col_ou <- "#0072b2"
col4 <- "#56b4e9"
col5 <- "#E6E6FA"


### K2 plot

dn <- df[c('SpName','V5','V6','new_names','Region')]
gk2 <- dn %>% gather(key = 'Cluster', value = 'Admixture', c(2,3))
gk2$SpName <- factor(gk2$SpName, levels = as.vector(df$SpName))

p2 <- ggplot(gk2) + aes(x = SpName, y = Admixture, fill = Cluster) +
  geom_bar(stat = 'identity', col=NA, width=1) +
  scale_fill_manual(values = c(col_ou,col4)) +
  scale_y_continuous(breaks = c(0,0.5,1), labels = c("0.00","0.50","1.00")) +
  labs(y = "K = 2") +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none",
        #axis.text.x = element_text(angle=90),
        axis.text.y = element_text(size=14),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_blank())
p2



### K3 plot

dn3 <- k3pop[c('SpName','V5','V6','V7','Region')]
gk3 <- dn3 %>% gather(key = 'Cluster', value = 'Admixture', c(2,3,4))
gk3$SpName <- factor(gk3$SpName, levels = as.vector(df$SpName))

p3 <- ggplot(gk3) + aes(x = SpName, y = Admixture, fill = Cluster) +
  geom_bar(stat = 'identity', col=NA, width=1) +
  scale_fill_manual(values = c(col_ou,col4,col5)) +
  scale_y_continuous(breaks = c(0,0.5,1), labels = c("0.00","0.50","1.00")) +
  labs(y = "K = 3") +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none",
        axis.text.y = element_text(size=14),
        axis.title.y = element_text(size=18),
        axis.text.x = element_blank())
p3




### Plotting regions

df$Stat <- 1.00
df$SpName <- factor(df$SpName, levels = df$SpName)

pr <- ggplot(df) + aes(x = SpName, y = Stat, fill = as.factor(Region)) +
  geom_bar(stat = 'identity', col=NA, width=1) +
  geom_text(x = 6, y = 0.5, label = "Europe", size = 7) + 
  geom_text(x = 16.5, y = 0.5, label = "North America", size = 7) + 
  scale_fill_manual(values = c("grey60","grey80")) +
  scale_y_continuous(breaks = c(0,1), labels = c("0.00","1.00")) +
  labs(y = "R") +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=18, colour = "white"),
        legend.position = "none",
        #axis.text.x = element_text(angle=90),
        #axis.text.x = element_text(size=14, angle=90),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=14, colour = "white"))
pr



# combine plots

grid.arrange(p2,p3,pr,ncol=1)

pdf(file="Figure_ulmi_structure_regions.pdf",width=6,height=2.5)
#lay <- rbind(c(1),c(1),c(1),c(2),c(2),c(2),c(3))
lay <- rbind(c(1),c(1),c(2),c(2),c(2),c(3),c(3),c(3))
grid.arrange(pr,p2,p3, ncol = 1, layout_matrix = lay)
dev.off()
