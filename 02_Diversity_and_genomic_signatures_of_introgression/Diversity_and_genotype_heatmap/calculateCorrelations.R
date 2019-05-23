library(dplyr)
library(ggplot2)
library(gridExtra)

### Reading file with stats

d <- read.table(file = "Windows_Statistics.txt", sep = "\t", header = TRUE)
stat = d$Dxy_A1U


####### PLOTS OF DXY ULM VS AME1, AME2 AND NOV

mytheme = theme(panel.background = element_blank(),
                axis.line = element_line(),
                axis.text = element_text(size=8),
                axis.title = element_text(size = 10))

### 1

rho = cor(d$Dxy_A1U,d$Dxy_A2U,method = "spearman")
sr = cor.test(d$Dxy_A1U,d$Dxy_A2U,method = "spearman")
pval = sr$p.value

p1 <- ggplot(d) + aes(x = Dxy_A1U, y = Dxy_A2U) +
  geom_abline(colour = "grey") +
  geom_point(fill = "#0072B2",colour = "black", shape = 21) +
  xlim(0.01,0.05) + ylim(0.01,0.05) +
  labs(x = "Dxy AME1 vs. ULM", y = "Dxy AME2 vs. ULM") +
  geom_text(x = 0.015, y = 0.045, label = paste("rho = ",round(rho,2)), hjust = 0.015, size = 4) +
  geom_text(x = 0.015, y = 0.041, label = paste("P-value < 2.2e-16"), hjust = 0.015, size = 4) +
  mytheme
p1

### 2

rho = cor(d$Dxy_A1U,d$Dxy_NU,method = "spearman")
sr = cor.test(d$Dxy_A1U,d$Dxy_NU,method = "spearman")
pval = sr$p.value

p2 <- ggplot(d) + aes(x = Dxy_A1U, y = Dxy_NU) +
  geom_abline(colour = "grey") +
  geom_point(fill = "#0072B2",colour = "black", shape = 21) +
  xlim(0.01,0.05) + ylim(0.01,0.05) +
  labs(x = "Dxy AME1 vs. ULM", y = "Dxy NOV vs. ULM") +
  geom_text(x = 0.015, y = 0.045, label = paste("rho = ",round(rho,2)), hjust = 0.015, size = 4) +
  geom_text(x = 0.015, y = 0.041, label = paste("P-value < 2.2e-16"), hjust = 0.015, size = 4) +
  mytheme
p2


### 3

rho = cor(d$Dxy_A2U,d$Dxy_NU,method = "spearman")
sr = cor.test(d$Dxy_A2U,d$Dxy_NU,method = "spearman")
pval = sr$p.value

p3 <- ggplot(d) + aes(x = Dxy_A2U, y = Dxy_NU) +
  geom_abline(colour = "grey") +
  geom_point(fill = "#0072B2",colour = "black", shape = 21) +
  xlim(0.01,0.05) + ylim(0.01,0.05) +
  labs(x = "Dxy AME2 vs. ULM", y = "Dxy NOV vs. ULM") +
  geom_text(x = 0.015, y = 0.045, label = paste("rho = ",round(rho,2)), hjust = 0.015, size = 4) +
  geom_text(x = 0.015, y = 0.041, label = paste("P-value < 2.2e-16"), hjust = 0.015, size = 4) +
  mytheme
p3


##### PLOTTING FST BETWEEN ULM AND OTHERS


### 1

rho = cor(d$Fst_A1U,d$Fst_A2U,method = "spearman")
sr = cor.test(d$Fst_A1U,d$Fst_A2U,method = "spearman")
pval = sr$p.value

p4 <- ggplot(d) + aes(x = Fst_A1U, y = Fst_A2U) +
  geom_abline(colour = "grey") +
  geom_point(fill = "#0072B2",colour = "black", shape = 21) +
  xlim(0.5,1) + ylim(0.5,1) +
  labs(x = "FST AME1 vs. ULM", y = "FST AME2 vs. ULM") +
  geom_text(x = 0.8, y = 0.65, label = paste("rho = ",round(rho,2)), hjust = 0.8, size = 4) +
  geom_text(x = 0.8, y = 0.61, label = paste("P-value < 2.2e-16"), hjust = 0.8, size = 4) +
  mytheme
p4

### 2

rho = cor(d$Fst_A1U,d$Fst_NU,method = "spearman")
sr = cor.test(d$Fst_A1U,d$Fst_NU,method = "spearman")
pval = sr$p.value

p5 <- ggplot(d) + aes(x = Fst_A1U, y = Fst_NU) +
  geom_abline(colour = "grey") +
  geom_point(fill = "#0072B2",colour = "black", shape = 21) +
  xlim(0.5,1) + ylim(0.5,1) +
  labs(x = "FST AME1 vs. ULM", y = "FST NOV vs. ULM") +
  geom_text(x = 0.8, y = 0.65, label = paste("rho = ",round(rho,2)), hjust = 0.8, size = 4) +
  geom_text(x = 0.8, y = 0.61, label = paste("P-value < 2.2e-16"), hjust = 0.8, size = 4) +
  mytheme
p5


### 3

rho = cor(d$Fst_A2U,d$Fst_NU,method = "spearman")
sr = cor.test(d$Fst_A2U,d$Fst_NU,method = "spearman")
pval = sr$p.value

p6 <- ggplot(d) + aes(x = Fst_A2U, y = Fst_NU) +
  geom_abline(colour = "grey") +
  geom_point(fill = "#0072B2",colour = "black", shape = 21) +
  xlim(0.5,1) + ylim(0.5,1) +
  labs(x = "FST AME2 vs. ULM", y = "FST NOV vs. ULM") +
  geom_text(x = 0.8, y = 0.65, label = paste("rho = ",round(rho,2)), hjust = 0.8, size = 4) +
  geom_text(x = 0.8, y = 0.61, label = paste("P-value < 2.2e-16"), hjust = 0.8, size = 4) +
  mytheme
p6


###### PLOTTING DXY VS PI

### 1

rho = cor(d$Dxy_A1U,d$Pi_AME1,method = "spearman")
sr = cor.test(d$Dxy_A1U,d$Pi_AME1,method = "spearman")
pval = sr$p.value

p7 <- ggplot(d) + aes(x = Dxy_A1U, y = Pi_AME1) +
  geom_point(fill = "#0072B2",colour = "black", shape = 21) +
  xlim(0,0.06) + ylim(0,0.02) +
  labs(x = "Dxy AME1 vs. ULM", y = "Nucl. diversity AME1") +
  geom_text(x = 0.036, y = 0.018, label = paste("rho = ",round(rho,2)), size = 4) +
  geom_text(x = 0.036, y = 0.016, label = paste("P-value =",round(pval,3)), size = 4) +
  mytheme
p7



### 2

rho = cor(d$Dxy_A2U,d$Pi_AME2,method = "spearman")
sr = cor.test(d$Dxy_A2U,d$Pi_AME2,method = "spearman")
pval = sr$p.value

p8 <- ggplot(d) + aes(x = Dxy_A2U, y = Pi_AME2) +
  geom_point(fill = "#0072B2",colour = "black", shape = 21) +
  xlim(0,0.06) + ylim(0,0.02) +
  labs(x = "Dxy AME2 vs. ULM", y = "Nucl. diversity AME2") +
  geom_text(x = 0.036, y = 0.018, label = paste("rho = ",round(rho,2)), size = 4) +
  geom_text(x = 0.036, y = 0.016, label = paste("P-value =",round(pval,3)), size = 4) +
  mytheme
p8


### 3

rho = cor(d$Dxy_NU,d$Pi_NOV,method = "spearman")
sr = cor.test(d$Dxy_NU,d$Pi_NOV,method = "spearman")
pval = sr$p.value

p9 <- ggplot(d) + aes(x = Dxy_NU, y = Pi_NOV) +
  geom_point(fill = "#0072B2",colour = "black", shape = 21) +
  xlim(0,0.06) + ylim(0,0.02) +
  labs(x = "Dxy NOV vs. ULM", y = "Nucl. diversity NOV") +
  geom_text(x = 0.036, y = 0.018, label = paste("rho = ",round(rho,2)), size = 4) +
  geom_text(x = 0.036, y = 0.016, label = paste("P-value =",round(pval,6)), size = 4) +
  mytheme
p9


######## PLOTTING  PI vs DXY intra-O. novo-ulmi

### ame1 - ame1/ame2
rho = cor(d$Dxy_A1A2,d$Pi_AME1,method = "spearman")
sr = cor.test(d$Dxy_A1A2,d$Pi_AME1,method = "spearman")
pval = sr$p.value

p10 <- ggplot(d) + aes(x = Dxy_A1A2, y = Pi_AME1) +
  geom_point(fill = "#0072B2",colour = "black", shape = 21) +
  xlim(0,0.06) + ylim(0,0.02) +
  labs(x = "Dxy AME1 vs. AME2", y = "Nucl. diversity AME1") +
  geom_text(x = 0.036, y = 0.018, label = paste("rho = ",round(rho,2)), size = 4) +
  geom_text(x = 0.036, y = 0.016, label = paste("P-value < 2.2e-16"), size = 4) +
  mytheme
p10


### ame1 - ame1/nov
rho = cor(d$Dxy_A1N,d$Pi_AME1,method = "spearman")
sr = cor.test(d$Dxy_A1N,d$Pi_AME1,method = "spearman")
pval = sr$p.value

p11 <- ggplot(d) + aes(x = Dxy_A1N, y = Pi_AME1) +
  geom_point(fill = "#0072B2",colour = "black", shape = 21) +
  xlim(0,0.06) + ylim(0,0.02) +
  labs(x = "Dxy AME1 vs. NOV", y = "Nucl. diversity AME1") +
  geom_text(x = 0.036, y = 0.018, label = paste("rho = ",round(rho,2)), size = 4) +
  geom_text(x = 0.036, y = 0.016, label = paste("P-value < 2.2e-16"), size = 4) +
  mytheme
p11


### ame2 - ame1/ame2
rho = cor(d$Dxy_A1A2,d$Pi_AME2,method = "spearman")
sr = cor.test(d$Dxy_A1A2,d$Pi_AME2,method = "spearman")
pval = sr$p.value

p12 <- ggplot(d) + aes(x = Dxy_A1A2, y = Pi_AME1) +
  geom_point(fill = "#0072B2",colour = "black", shape = 21) +
  xlim(0,0.06) + ylim(0,0.02) +
  labs(x = "Dxy AME1 vs. AME2", y = "Nucl. diversity AME2") +
  geom_text(x = 0.036, y = 0.018, label = paste("rho = ",round(rho,2)), size = 4) +
  geom_text(x = 0.036, y = 0.016, label = paste("P-value < 2.2e-16"), size = 4) +
  mytheme
p12

### ame2 - ame2/nov
rho = cor(d$Dxy_A2N,d$Pi_AME2,method = "spearman")
sr = cor.test(d$Dxy_A2N,d$Pi_AME2,method = "spearman")
pval = sr$p.value

p13 <- ggplot(d) + aes(x = Dxy_A2N, y = Pi_AME2) +
  geom_point(fill = "#0072B2",colour = "black", shape = 21) +
  xlim(0,0.06) + ylim(0,0.02) +
  labs(x = "Dxy AME2 vs. NOV", y = "Nucl. diversity AME2") +
  geom_text(x = 0.036, y = 0.018, label = paste("rho = ",round(rho,2)), size = 4) +
  geom_text(x = 0.036, y = 0.016, label = paste("P-value < 2.2e-16"), size = 4) +
  mytheme
p13


### nov - ame1/nov
rho = cor(d$Dxy_A1N,d$Pi_NOV,method = "spearman")
sr = cor.test(d$Dxy_A1N,d$Pi_NOV,method = "spearman")
pval = sr$p.value

p14 <- ggplot(d) + aes(x = Dxy_A1N, y = Pi_NOV) +
  geom_point(fill = "#0072B2",colour = "black", shape = 21) +
  xlim(0,0.06) + ylim(0,0.02) +
  labs(x = "Dxy AME1 vs. NOV", y = "Nucl. diversity NOV") +
  geom_text(x = 0.036, y = 0.018, label = paste("rho = ",round(rho,2)), size = 4) +
  geom_text(x = 0.036, y = 0.016, label = paste("P-value < 2.2e-16"), size = 4) +
  mytheme
p14


### nov- ame2/nov
rho = cor(d$Dxy_A2N,d$Pi_NOV,method = "spearman")
sr = cor.test(d$Dxy_A2N,d$Pi_NOV,method = "spearman")
pval = sr$p.value

p15 <- ggplot(d) + aes(x = Dxy_A2N, y = Pi_NOV) +
  geom_point(fill = "#0072B2",colour = "black", shape = 21) +
  xlim(0,0.06) + ylim(0,0.02) +
  labs(x = "Dxy AME2 vs. NOV", y = "Nucl. diversity NOV") +
  geom_text(x = 0.036, y = 0.018, label = paste("rho = ",round(rho,2)), size = 4) +
  geom_text(x = 0.036, y = 0.016, label = paste("P-value < 2.2e-16"), size = 4) +
  mytheme
p15


### Putting all plots together

pdf('calculateCorrelations.pdf',w=8,h=10)
grid.arrange(p1,p2,p3,
             p4,p5,p6,
             p7,p8,p9,
             p10,p12,p14,
             p11,p13,p15,
             ncol=3)
dev.off()
