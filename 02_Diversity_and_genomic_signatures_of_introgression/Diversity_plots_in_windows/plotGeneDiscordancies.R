rm(list=ls())
source("plot_twisst.R")
library(gridExtra)
library(dplyr)
library(ggplot2)
library(lemon)
library(grid)
library(tidyr)


###########################################################
############### PLOTTING TREE WEIGHTS #####################
###########################################################


### 100 SNPs window

weights_file <- "./data/combineTwisstOutput_topos.out"
window_data_file <- "./data/combineTwisstOutput_coord.out"
w = read.table(weights_file, header = TRUE)
weights <- w[,c(1,13,4,3,2,9,6,12,5,15,10,7,14,11,8)]
colnames(weights) <- c("topo1","topo2","topo3","topo4","topo5","topo6","topo7","topo8","topo9","topo10","topo11","topo12","topo13","topo14","topo15")
order <- c("topo1","topo2","topo3","topo4","topo5","topo6","topo7","topo8","topo9","topo10","topo11","topo12","topo13","topo14","topo15")
#weights <- weights / apply(weights, 1, sum)
topoNames = names(weights)
window_data = read.table(window_data_file, header = TRUE)
good_rows = which(is.na(apply(weights,1,sum)) == FALSE)
weights <- weights[good_rows,]
window_data = window_data[good_rows,]
cols <- c("#6699CC","#DDCC77","#CC6677","#44AA99","#999933","maroon1","#BF864D","#86BF4D","purple","red","limegreen","#4D4DBF","#BF4D86","#FFD39B","#CDC8B1")
rest="black"
cols4 <- c("#6699CC","#DDCC77","#CC6677",rest,rest,rest,rest,rest,rest,rest,rest,rest,rest,rest,rest)

weights_smooth = smooth.weights(window_positions=window_data$mid, weights_dataframe = weights,
                                span = 0.008, window_sites=window_data$sites)

ws <- weights_smooth[,c(15:1)]
tot <- cbind(ws,window_data)
T <- tot %>% gather("Topo","Prop",topo15:topo1)
T$Topo <- factor(T$Topo, levels = order)

pT <- ggplot(T) +
  aes(x = mid, y = Prop, fill = Topo) +
  geom_area() +
  scale_fill_manual(values = cols4) +
  #labs(y = "T. w.")+
  ggtitle("Tree weighting 100 SNPs") +
  theme_bw() +
  theme(axis.title=element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(),
        axis.line.y = element_line(),
        legend.position="none") +
  geom_vline(xintercept = c(6937932,13755643,17425415,20845118,23693821,26495415,29253639,31784886), 
             colour = "white", size=1)+
  scale_x_continuous(limits = c(0, 32000000),
                     breaks=c(3468966,10346787.5,15590529,19135266.5,22269469.5,25094618,27874527,30519262.5),
                     labels = c("1","2","3","4","5","6","7","8")) +
  coord_capped_cart(xlim = c(0, 32000000),bottom='both',left='both')
#pT


### 100 SNPs windows with >=3 SNPs in OU and ONU

weights_file <- "./data/combineTwisstOutput_min3_topos.out"
window_data_file <- "./data/combineTwisstOutput_min3_coord.out"
w = read.table(weights_file, header = TRUE)
weights <- w[,c(1,13,4,3,2,9,6,12,5,15,10,7,14,11,8)]
colnames(weights) <- c("topo1","topo2","topo3","topo4","topo5","topo6","topo7","topo8","topo9","topo10","topo11","topo12","topo13","topo14","topo15")
order <- c("topo1","topo2","topo3","topo4","topo5","topo6","topo7","topo8","topo9","topo10","topo11","topo12","topo13","topo14","topo15")
#weights <- weights / apply(weights, 1, sum)
topoNames = names(weights)
window_data = read.table(window_data_file, header = TRUE)
good_rows = which(is.na(apply(weights,1,sum)) == FALSE)
weights <- weights[good_rows,]
window_data = window_data[good_rows,]
cols <- c("#6699CC","#DDCC77","#CC6677","#44AA99","#999933","maroon1","#BF864D","#86BF4D","purple","red","limegreen","#4D4DBF","#BF4D86","#FFD39B","#CDC8B1")
rest="black"
cols4 <- c("#6699CC","#DDCC77","#CC6677",rest,rest,rest,rest,rest,rest,rest,rest,rest,rest,rest,rest)

weights_smooth = smooth.weights(window_positions=window_data$mid, weights_dataframe = weights,
                                span = 0.008, window_sites=window_data$sites)

ws <- weights_smooth[,c(15:1)]
tot <- cbind(ws,window_data)
T <- tot %>% gather("Topo","Prop",topo15:topo1)
T$Topo <- factor(T$Topo, levels = order)

pT2 <- ggplot(T) +
  aes(x = mid, y = Prop, fill = Topo) +
  geom_area() +
  scale_fill_manual(values = cols4) +
  #labs(y = "T. w.")+
  ggtitle("Tree weighting 100 SNPs with >= 3 SNPs in OU and ONU") +
  theme_bw() +
  theme(axis.title=element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(),
        axis.line.y = element_line(),
        legend.position="none") +
  geom_vline(xintercept = c(6937932,13755643,17425415,20845118,23693821,26495415,29253639,31784886), 
             colour = "white", size=1)+
  scale_x_continuous(limits = c(0, 32000000),
                     breaks=c(3468966,10346787.5,15590529,19135266.5,22269469.5,25094618,27874527,30519262.5),
                     labels = c("1","2","3","4","5","6","7","8")) +
  coord_capped_cart(xlim = c(0, 32000000),bottom='both',left='both')
#pT2





### 50 kb windows

weights_file <- "./data/combineTwisstOutput_50k_topos_head.out"
window_data_file <- "./data/combineTwisstOutput_50k_coord.out"
w = read.table(weights_file, header = TRUE)
weights <- w[,c(1,13,4,3,2,9,6,12,5,15,10,7,14,11,8)]
colnames(weights) <- c("topo1","topo2","topo3","topo4","topo5","topo6","topo7","topo8","topo9","topo10","topo11","topo12","topo13","topo14","topo15")
order <- c("topo1","topo2","topo3","topo4","topo5","topo6","topo7","topo8","topo9","topo10","topo11","topo12","topo13","topo14","topo15")
weights <- weights / apply(weights, 1, sum)
topoNames = names(weights)
window_data = read.table(window_data_file, header = TRUE, sep = "\t")
good_rows = which(is.na(apply(weights,1,sum)) == FALSE)
weights <- weights[good_rows,]
window_data = window_data[good_rows,]
cols <- c("#6699CC","#DDCC77","#CC6677","#44AA99","#999933","maroon1","#BF864D","#86BF4D","purple","red","limegreen","#4D4DBF","#BF4D86","#FFD39B","#CDC8B1")
rest="black"
cols4 <- c("#6699CC","#DDCC77","#CC6677",rest,rest,rest,rest,rest,rest,rest,rest,rest,rest,rest,rest)


weights_smooth = smooth.weights(window_positions=window_data$mid, weights_dataframe = weights,
                                span = 0.006, window_sites=window_data$sites)

ws <- weights_smooth[,c(15:1)]
tot <- cbind(ws,window_data)
T <- tot %>% gather("Topo","Prop",topo15:topo1)
T$Topo <- factor(T$Topo, levels = order)

pT3 <- ggplot(T) +
  aes(x = mid, y = Prop, fill = Topo) +
  geom_area() +
  scale_fill_manual(values = cols4) +
  #labs(y = "T. w.")+
  ggtitle("Tree weighting 50kb") +
  theme_bw() +
  theme(axis.title=element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(),
        axis.line.y = element_line(),
        legend.position="none") +
  geom_vline(xintercept = c(6937932,13755643,17425415,20845118,23693821,26495415,29253639,31784886), 
             colour = "white", size=1)+
  scale_x_continuous(limits = c(0, 32000000),
                     breaks=c(3468966,10346787.5,15590529,19135266.5,22269469.5,25094618,27874527,30519262.5),
                     labels = c("1","2","3","4","5","6","7","8")) +
  coord_capped_cart(xlim = c(0, 32000000),bottom='both',left='both')
#pT3






### Plotting all graphs together

pdf(file = paste("Figure_TreeWeighting_different_windows.pdf",sep=""),width=8.27,height=4)
grid.arrange(arrangeGrob(pT,
                         pT2,
                         pT3,
                         ncol=1,
                         bottom = textGrob("Chromosome", hjust = 0.1)))
dev.off()







