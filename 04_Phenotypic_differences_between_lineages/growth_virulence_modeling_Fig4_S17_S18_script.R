library(emmeans)
library(ggplot2)
library(ggsignif)

##### Growth phenotypes #####
d=read.csv2("data_growth_experiment_introgression_infos.csv")
str(d)

### Best model (data exporation not shown) ###

bestmodel=lmer(sqrt(slope) ~ lineage + media + temperature + PC1 + PC3 + lineage:temperature + lineage:media + media:temperature + temperature:PC1 + (1|rep), data = d)
residplot(bestmodel)
Anova(bestmodel)
summary(bestmodel)

### Model followed by post hoc tests ###
fm=lmer(sqrt(slope)~media+lineage+temperature+lineage:temperature+lineage:media+media:temperature+(1|rep), data=d)
emmeans(fm,list(pairwise~lineage|media+temperature),adjust="tukey")

### Effect of introgression ###
d_ONU=subset(d, lineage!="OU")
d_ONU=droplevels(d_ONU)

fm=lmer(sqrt(slope) ~ lineage + temperature + media + int + chr1 + chr2 + chr6 + chr8 + lineage:temperature + lineage:media + lineage:int + temperature:media + temperature:int + media:int + (1|rep), data=d_ONU)
residplot(fm)
Anova(fm)
summary(fm)
emmeans(fm,list(pairwise~int|temperature),adjust="tukey")

##### Virulence on apples #####
apple=read.table("data_apple_int.txt", head=T)
apple_ONU=subset(apple, lineage != "ulm")

lm=lmer(pred~as.factor(chr1)+as.factor(chr2)+as.factor(chr6)+as.factor(chr8)+int+lineage+(1|pomme_num), data=apple_ONU)
residplot(lm)
Anova(lm)
summary(lm)

##### Plots #####

### Fig 4a ###
ggplot(d, aes(PC1,slope/1000, color=temperature)) +
  geom_point() + geom_smooth(method="lm") +
  scale_color_manual(values=c("yellow","orange","red")) +
  theme_bw() +
  scale_y_continuous(labels= scales::comma, limits=c(0,1200)) +
  ylab("Growth rate (in 1000 pixels/hour)") +
  xlab("PC1") +
  theme(axis.text = element_text(size=28), 
        axis.title=element_text(size=30,face="bold"), 
        axis.title.y = element_text(margin = margin(t=20,r=20,b=20,l=20)), 
        axis.title.x = element_text(margin = margin(t=20,r=20,b=20,l=20)),
        legend.text=element_text(size=30),
        legend.title = element_text(size=30),
        plot.margin = unit(c(1,1,1.5,1.2),"cm"))

ggplot(d, aes(PC2,slope/1000)) +
  geom_point() + geom_smooth(method="lm") +
  theme_bw() +
  scale_y_continuous(labels= scales::comma, limits=c(0,1200)) +
  ylab("Growth rate (in 1000 pixels/hour)") +
  xlab("PC2") +
  theme(axis.text = element_text(size=30), 
        axis.title=element_text(size=30,face="bold"), 
        axis.title.y = element_text(margin = margin(t=20,r=20,b=20,l=20)), 
        axis.title.x = element_text(margin = margin(t=20,r=20,b=20,l=20)),
        plot.margin = unit(c(1,1,1.5,1.2),"cm"))

ggplot(d, aes(PC3,slope/1000)) +
  geom_point() + geom_smooth(method="lm") +
  theme_bw() +
  scale_y_continuous(labels= scales::comma, limits=c(0,1200)) +
  ylab("Growth rate (in 1000 pixels/hour)") +
  xlab("PC3") +
  theme(axis.text = element_text(size=30), 
        axis.title=element_text(size=30,face="bold"), 
        axis.title.y = element_text(margin = margin(t=20,r=20,b=20,l=20)), 
        axis.title.x = element_text(margin = margin(t=20,r=20,b=20,l=20)),
        plot.margin = unit(c(1,1,1.5,1.2),"cm"))

### Fig 4b ###
d_ONU_MEA=subset(d_ONU, media=="MEA")

plot=ggplot(data=(d_ONU_MEA),aes(x = as.factor(int), y = as.numeric(as.character(slope/1000)))) + 
  geom_boxplot(aes(fill=as.factor(int)),lwd=0.3, width=0.65) + 
  facet_grid(temperature~.) + 
  theme_bw() + 
  scale_fill_manual(values=c("#f0d9d9","#cd8382")) + 
  ylab("Growth rate (in 1000 pixels/hour)") + xlab("Introgression from OU") +
  theme(axis.line = element_line(colour = 'black', size = 0.2), 
        strip.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        axis.text.y = element_text(colour="black"),
        axis.text.x = element_text(colour="black"),
        axis.text=element_text(size=18),
        axis.title=element_text(size=20,face="bold"),
        strip.text.x = element_text(size = 18), 
        strip.text.y = element_text(size = 18),
        axis.title.y = element_text(margin = margin(t=20,r=20,b=20,l=20)), 
        axis.title.x = element_text(margin = margin(t=20,r=20,b=20,l=20))) + 
  theme(legend.position="none", panel.spacing = unit(1,"lines"), plot.margin = unit(c(1,1,1.5,1.2),"cm")) +
  scale_y_continuous(labels= scales::comma, limits=c(0,1200))#+

plot + geom_signif(comparisons = list(c("0","1")), map_signif_level = T)


### Fig 4c ###
# Chromosome 1
plot_chr1=ggplot(data=na.omit(apple_ONU),aes(x = as.factor(chr1), y = as.numeric(as.character(pred)))) + 
  geom_boxplot(aes(fill=as.factor(chr1)),lwd=0.3) + 
  ylab("Necrosis size on apple (in mm)") + xlab("Introgression on chromosome 1") +
  scale_fill_manual(values=c("#f0d9d9","#cd8382")) +
  theme(axis.line = element_line(colour = 'black', size = 0.5), 
        strip.background = element_blank(), 
        panel.background = element_blank(),
        axis.text.y = element_text(colour="black"),
        axis.text.x = element_text(colour="black"),
        axis.text=element_text(size=26),
        axis.title=element_text(size=28,face="bold"),
        strip.text.x = element_text(size = 26), 
        strip.text.y = element_text(size = 26),
        axis.title.y = element_text(margin = margin(t=20,r=20,b=20,l=20)), 
        axis.title.x = element_text(margin = margin(t=20,r=20,b=20,l=20)))+
  scale_y_continuous(limits=c(10,36), breaks = seq(0,36,5))+
  theme(legend.position="none", plot.margin = unit(c(1.5,1.5,1.5,1.5),"cm")) #+

df1 <- data.frame(a = c(1,1,2,2,2), b = c(34.5, 35, 35, 35, 34.5))

plot_chr1 + geom_line(data = df1, aes(x = a, y = b), size=1) + annotate("text", x = 1.5, y = 35.5, label = "*", size = 18)

# Chromosome 8
plot_chr8=ggplot(data=na.omit(apple_ONU),aes(x = as.factor(chr8), y = as.numeric(as.character(pred)))) + 
  geom_boxplot(aes(fill=as.factor(chr8)),lwd=0.3) + 
  ylab("Necrosis size on apple (in mm)") + xlab("Introgression on chromosome 8") +
  scale_fill_manual(values=c("#f0d9d9","#cd8382")) +
  theme(axis.line = element_line(colour = 'black', size = 0.5), 
        strip.background = element_blank(), 
        panel.background = element_blank(),
        axis.text.y = element_text(colour="black"),
        axis.text.x = element_text(colour="black"),
        axis.text=element_text(size=26),
        axis.title=element_text(size=28,face="bold"),
        strip.text.x = element_text(size = 26), 
        strip.text.y = element_text(size = 26),
        axis.title.y = element_text(margin = margin(t=20,r=20,b=20,l=20)), 
        axis.title.x = element_text(margin = margin(t=20,r=20,b=20,l=20)))+
  scale_y_continuous(limits=c(10,36), breaks = seq(0,36,5))+
  theme(legend.position="none", plot.margin = unit(c(1.5,1.5,1.5,1.5),"cm")) 

df2 <- data.frame(a = c(1,1,2,2,2), b = c(34.5, 35, 35, 35, 34.5))

plot_chr8 + geom_line(data = df2, aes(x = a, y = b), size=1) + annotate("text", x = 1.5, y = 35.5, label = "***", size = 18)


### Fig S17 ###

d$group=as.factor(as.character(ifelse(d$lineage=="AME2"&d$int==1, "AME2_int",ifelse(d$lineage=="AME1"&d$int==1, "AME1_int",ifelse(d$lineage=="NOV"&d$int==1, "NOV_int",as.character(d$lineage))))))


plot=ggplot(data=na.omit(d),aes(x = as.factor(group), y = as.numeric(as.character(slope/1000)))) + 
  geom_boxplot(aes(fill=as.factor(group)),lwd=0.3) + 
  facet_grid(media~temperature) + 
  theme_bw() + 
  scale_fill_manual(values=c("#d55e00","#954100","#f0e442","#a89f2e","#009e73","#006e50","#0072b2")) + 
  ylab("Growth rate (in 1000 pixels/h)") + xlab("Lineage") +
  theme(axis.line = element_line(colour = 'black', size = 0.2), 
        strip.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        axis.text.y = element_text(colour="black"),
        axis.text.x = element_text(colour="black", angle=90, hjust=1),
        axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"),
        strip.text.x = element_text(size = 16), 
        strip.text.y = element_text(size = 16),
        axis.title.y = element_text(margin = margin(t=20,r=20,b=20,l=20)), 
        axis.title.x = element_text(margin = margin(t=20,r=20,b=20,l=20))) + 
  theme(legend.position="none") + scale_y_continuous(labels= scales::comma, limits=c(0,1200))

plot

### Fig S18 ###
apple$species=ifelse(apple$lineage=="ulm",as.character("OU"),as.character("ONU"))

plot=ggplot(data=na.omit(apple),aes(x = as.factor(species), y = pred)) + 
  geom_boxplot(aes(fill=as.factor(species)),lwd=0.3, width=0.65) + 
  theme_bw() + 
  scale_fill_manual(values=c("#DEDBB1","#0072b2")) + 
  ylab("Necrosis size (in mm)") + xlab("Species") +
  theme(axis.line = element_line(colour = 'black', size = 0.2), 
        strip.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        axis.text.y = element_text(colour="black"),
        axis.text.x = element_text(colour="black"),
        axis.text=element_text(size=18),
        axis.title=element_text(size=20,face="bold"),
        strip.text.x = element_text(size = 18), 
        strip.text.y = element_text(size = 18),
        axis.title.y = element_text(margin = margin(t=20,r=20,b=20,l=20)), 
        axis.title.x = element_text(margin = margin(t=20,r=20,b=20,l=20))) + 
  theme(legend.position="none", panel.spacing = unit(1,"lines"), plot.margin = unit(c(1,1,1.5,1.2),"cm")) +
  scale_y_continuous(labels= scales::comma, limits=c(10,36))#+

plot + geom_signif(comparisons = list(c("OU","ONU")), map_signif_level = T, textsize = 7)






