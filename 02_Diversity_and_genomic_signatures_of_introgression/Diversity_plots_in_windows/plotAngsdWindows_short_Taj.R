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
  ggtitle("Tree weighting") +
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


###########################################################
################      PLOTTING DXY     ####################
###########################################################


pairs <- c("ame1.ame2","nov.ame1","nov.ame2")
dd <- data.frame()
for (n in 1:length(pairs)) {
  sp1 <- strsplit(pairs[n],split='[.]')[[1]][1]
  sp2 <- strsplit(pairs[n],split='[.]')[[1]][2]
  df <- read.table(file=paste("./data/",sp1,".",sp2,"_circos_Dxy.txt",sep=""),sep="\t",header=FALSE)
  colnames(df) <- c("Chr","Midpoint","Midpoint2","Stat")
  df["Pair"] <- paste(sp1,".",sp2,sep="")
  dd <- rbind(dd,df)
}

P <- matrix(nrow = nrow(dd), ncol = 1)
for (i in 1:nrow(dd)) {
  r = dd[i,]
  if (r[1] == "OphioH327chr_1") {
    p = r[2] + 0
    P[i,1] = p[1,1]}
  else if (r[1] == "OphioH327chr_2") {
    p = r[2] + 6937932
    P[i,1] = p[1,1]}
  else if (r[1] == "OphioH327chr_3") {
    p = r[2] + 13755643
    P[i,1] = p[1,1]}
  else if (r[1] == "OphioH327chr_4") {
    p = r[2] + 17425415
    P[i,1] = p[1,1]}
  else if (r[1] == "OphioH327chr_5") {
    p = r[2] + 20845118
    P[i,1] = p[1,1]}
  else if (r[1] == "OphioH327chr_6") {
    p = r[2] + 23693821
    P[i,1] = p[1,1]}
  else if (r[1] == "OphioH327chr_7") {
    p = r[2] + 26495415
    P[i,1] = p[1,1]}
  else if (r[1] == "OphioH327chr_8") {
    p = r[2] + 29253639
    P[i,1] = p[1,1]}
}
dd["Position"] <- P[,1]


plot_dxy = list()
for (pair in 1:length(pairs)) {
  pop1 <- gsub("ULM","OU",toupper(unlist(strsplit(pairs[pair],"[.]"))[1]))
  pop2 <- gsub("ULM","OU",toupper(unlist(strsplit(pairs[pair],"[.]"))[2]))
  dd_pair <- dd %>% filter(Pair == pairs[pair])
  p <-  ggplot(dd_pair) +
    aes(x = Position, y = Stat, fill = Pair, colour = Chr) +
    geom_point(size=1) +
    scale_colour_manual(values = c("black","grey65","black","grey65","black","grey65","black","grey65")) +
    scale_x_continuous(limits = c(0, 32000000),
                       breaks=c(3468966,10346787.5,15590529,19135266.5,22269469.5,25094618,27874527,30519262.5),
                       labels = c("1","2","3","4","5","6","7","8")) +
    scale_y_continuous(limits = c(0,0.05),breaks = c(0.0,0.025,0.05),
                       labels = c("0.000","0.025","0.050")) +
    ggtitle(paste("Dxy",pop1,"vs.",pop2,sep=" ")) +
    theme_bw() +
    theme(axis.title=element_blank(),
          plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line.x = element_line(),
          axis.line.y = element_line(),
          legend.position="None") +
    coord_capped_cart(xlim = c(0, 32000000),bottom='both',left='both')
  plot_dxy[[pair]] = p
}
plot_dxy[[1]]


###########################################################
################       PLOTTING Tajima     ################
###########################################################


species <- c("ame1","ame2","nov","ulm")
dp <- data.frame()
for (n in 1:length(species)) {
  df <- read.table(file=paste("./data/",species[n],"_circos_Tajima.txt",sep=""),sep="\t",header=FALSE)
  colnames(df) <- c("Chr","Midpoint","Midpoint2","Stat")
  df["Species"] <- species[n]
  dp <- rbind(dp,df)
}

R <- matrix(nrow = nrow(dp), ncol = 1)
for (i in 1:nrow(dp)) {
  r = dp[i,]
  if (r[1] == "OphioH327chr_1") {
    p = r[2] + 0
    R[i,1] = p[1,1]}
  else if (r[1] == "OphioH327chr_2") {
    p = r[2] + 6937932
    R[i,1] = p[1,1]}
  else if (r[1] == "OphioH327chr_3") {
    p = r[2] + 13755643
    R[i,1] = p[1,1]}
  else if (r[1] == "OphioH327chr_4") {
    p = r[2] + 17425415
    R[i,1] = p[1,1]}
  else if (r[1] == "OphioH327chr_5") {
    p = r[2] + 20845118
    R[i,1] = p[1,1]}
  else if (r[1] == "OphioH327chr_6") {
    p = r[2] + 23693821
    R[i,1] = p[1,1]}
  else if (r[1] == "OphioH327chr_7") {
    p = r[2] + 26495415
    R[i,1] = p[1,1]}
  else if (r[1] == "OphioH327chr_8") {
    p = r[2] + 29253639
    R[i,1] = p[1,1]}
}
dp["Position"] <- R[,1]



plot_td = list()
for (sp in 1:length(species)) {
  SP <- gsub("ULM","OU",toupper(species[sp]))
  dp_sp <- dp %>% filter(Species == species[sp])
  p <-  ggplot(dp_sp) +
    aes(x = Position, y = Stat, fill = Species, colour = Chr) +
    geom_hline(yintercept = 0, colour = "lightgrey") +
    geom_point(size=1) +
    scale_colour_manual(values = c("black","grey65","black","grey65","black","grey65","black","grey65")) +
    scale_x_continuous(limits = c(0, 32000000),
                       breaks=c(3468966,10346787.5,15590529,19135266.5,22269469.5,25094618,27874527,30519262.5),
                       labels = c("1","2","3","4","5","6","7","8")) +
    scale_y_continuous(limits = c(-4.0,4.0),breaks = c(-4.0,0.0,4.0),
                       labels = c("-4.00","0","4.00")) +
    ggtitle(paste("Tajima's D",SP,sep=" ")) +
    theme_bw() +
    theme(axis.title=element_blank(),
          plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line.x = element_line(),
          axis.line.y = element_line(),
          legend.position="None") +
    coord_capped_cart(xlim = c(0, 32000000),bottom='both',left='both')
  plot_td[[sp]] = p
}

plot_td[[1]]


### Plotting all graphs together


pdf(file = paste("Figure_TreeWeighting_Dxy_Tajima.pdf",sep=""),width=8.27,height=8)
grid.arrange(arrangeGrob(pT,
                         plot_dxy[[1]],plot_dxy[[2]],plot_dxy[[3]],
                         plot_td[[1]],plot_td[[2]],plot_td[[3]],plot_td[[4]],
                         ncol=1,
                         bottom = textGrob("Chromosome", hjust = 0.1)))
dev.off()


