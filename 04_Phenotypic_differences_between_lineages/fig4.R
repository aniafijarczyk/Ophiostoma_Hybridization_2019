d=read.csv2("data_fig4.csv")

library(ggplot2)

plot=ggplot(data=na.omit(d),aes(x = lineage, y = slope)) + 
  geom_boxplot(aes(fill=lineage),lwd=0.3) + 
  facet_grid(media~temperature) + 
  theme_bw() + 
  scale_fill_manual(values=c("#CD2626","#EEB422","#698B22","#1874CD")) + 
  ylab("Growth rate (in pixel intensity)") + xlab("Lineage") +
  theme(axis.line = element_line(colour = 'black', size = 0.2), 
        strip.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        axis.text.y = element_text(colour="black"),
        axis.text.x = element_text(colour="black"),
        axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"),
        strip.text.x = element_text(size = 16), 
        strip.text.y = element_text(size = 16),
        axis.title.y = element_text(margin = margin(t=20,r=20,b=20,l=20)), 
        axis.title.x = element_text(margin = margin(t=20,r=20,b=20,l=20))) + 
  theme(legend.position="none") + scale_y_continuous(labels= scales::comma, limits=c(0,1100000))

plot
