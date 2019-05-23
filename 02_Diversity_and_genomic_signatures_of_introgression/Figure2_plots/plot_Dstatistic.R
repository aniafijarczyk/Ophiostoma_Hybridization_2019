# Script for plotting D-statistics

library(ggplot2)


### Plotting D-statistics

m <- data.frame(rbind(c("D(out,ulm,ame,nov)","D(H,U;A,N)"),c("D(out,ulm,ame1,nov)","D(H,U;A1,N)"),
           c("D(out,ulm,ame2,nov)","D(H,U;A2,N)"), c("D(out,nov,ame1,ame2)","D(H,N;A1,A2)"),
           c("D(ulm,nov,ame1,ame2)","D(U,N;A1,A2)"), c("D(out,ulm,ame1,ame2)","D(H,U;A1,A2)")))
colnames(m) <- c("Test","Short")

df <- read.table("Dstats_results.txt",sep="\t",header=TRUE)

df$Test <- factor(df$Test, levels = c("D(out,ulm,ame1H,ame1)","D(out,ulm,ame2H,ame2)","D(out,ulm,novH,nov)",
                                      "D(out,ame2,ame1H,ame1)","D(out,nov,ame1H,ame1)","D(out,ame1,ame2H,ame2)",
                                      "D(out,nov,ame2H,ame2)","D(out,ame1,novH,nov)","D(out,ame2,novH,nov)"))

p <- ggplot(df) + aes(x = Test, y = d) +
  geom_hline(yintercept = 0, colour = "grey") +
  geom_point(size=2) +
  geom_errorbar(aes(ymin = d - 3*se, ymax = d + 3*se),width = 0.3) +
  labs(y = expression(paste(italic(D),"-statistic")))+
  scale_x_discrete(labels = c("D(out,ulm,ame1H,ame1)" = expression(paste("ULM",'->',"AME1")),
                              "D(out,ulm,ame2H,ame2)" = expression(paste("ULM",'->',"AME2")),
                              "D(out,ulm,novH,nov)" = expression(paste("ULM",'->',"NOV")),
                              "D(out,ame2,ame1H,ame1)" = expression(paste("AME2",'->',"AME1")),
                              "D(out,nov,ame1H,ame1)" = expression(paste("NOV",'->',"AME1")),
                              "D(out,ame1,ame2H,ame2)" = expression(paste("AME1",'->',"AME2")),
                              "D(out,nov,ame2H,ame2)" = expression(paste("NOV",'->',"AME2")),
                              "D(out,ame1,novH,nov)" = expression(paste("AME1",'->',"NOV")),
                              "D(out,ame2,novH,nov)" = expression(paste("AME2",'->',"NOV")))) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size=16,angle=90),
        axis.text.y = element_text(size=16),
        axis.title.y = element_text(size=20),
        axis.title.x = element_blank())
p




pdf("plot_Dstatistic.pdf",w=10,h=6)
p
dev.off()


