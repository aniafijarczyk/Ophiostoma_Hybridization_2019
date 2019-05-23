
library(dplyr)
library(ggplot2)
library(ggforce)


### Plotting general diversity distribution (across windows) for each lineage

d <- read.table(file="./Windows_Statistics.txt",sep="\t",header=TRUE)
dt <- d %>% gather(key=SP,value=Pi,c(Pi_ULM,Pi_NOV,Pi_AME1,Pi_AME2))

ggplot(dt) + aes(x = SP, y = Pi) +
  geom_sina(size=0.5) +
  geom_boxplot(width=0.1,outlier.shape = NA) +
  scale_x_discrete(breaks = c("Pi_AME1","Pi_AME2","Pi_NOV","Pi_ULM"), labels = c("AME1","AME2","NOV","ULM")) +
  labs(x = "Lineage", y = "Nucleotide diversity") +
  theme_bw() + 
  theme(panel.grid = element_blank(),
        axis.text = element_text(size=18),
        axis.title = element_text(size=22))

