library(dplyr)
library(ggplot2)



### getting table with snipre results
d <- read.table("Input_snipre_8379_bayesian_results_50k.csv",sep=",",header=TRUE)
cols <- c( 'BSnIPRE.class' , 'BSnIPRE.Rclass' )
d$Class <- apply(d[ , cols ] , 1 , paste , collapse = "_" )


### Plotting results
p <- ggplot(d) + aes(x = BSnIPRE.Rest, y = BSnIPRE.est, colour = Class) +
  geom_hline(yintercept = 0, linetype = 2,colour = "grey") +
  geom_vline(xintercept = 0, linetype = 2,colour = "grey") +
  geom_point(size=2) +
  scale_colour_manual(values = c("#0072B2","darkgrey","#D55E00","#F0E442"), labels = c("Constraint","Neutral","Selection & constraint","Selection")) +
  labs(x = "Constraint effect", y = "Selection effect") +
  theme(panel.background = element_blank(),
        axis.text = element_text(size=10),
        axis.title = element_text(size=14),
        panel.border = element_rect(fill = NA),
        legend.key = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=12),
        legend.position = "top")
p


pdf('plotSnipreResults.pdf',w=10,h=6)
p
dev.off()
