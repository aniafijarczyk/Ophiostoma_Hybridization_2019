"
Script reads a table with SNPs and genotypes for multiple samples,
plots a heatmap with genotypes and performs clustering with hclust.
It also reads a file with sample names in the same order as they appear in the table.

USAGE:

Rscript plotGenotypesHeatMap.R

requirements:   convert2MajorMinor.tab
                strains_vcf_order.txt
output:         plotGenotypesHeatmap.pdf
                strains_hclust_order.txt
"



library(dplyr)
library(gplots)


# Reading files with genotypes and names

d <- read.table("convert2MajorMinor.tab",sep="\t",header=FALSE)
names <- read.table("strains_vcf_order.txt", sep="\t",header = TRUE)
namesv <- as.vector(names$SamName)
colnames(d) <- c("chrom","pos","num",namesv)


# Converting dataframe to a matrix

df <- d[seq(1, nrow(d), 300), ]
sub <- df[,c(4:97)]
rnames <- df[,c(3)]
data <- data.matrix(sub)
rownames(data) <- rnames
mydata <- t(data)


# Getting pop names and colors

popnames <- rownames(mydata)
popsymbols = c()
colsymbols = c()
for (i in 1:length(popnames)) {
  fields = strsplit(popnames[i],"_")
  pop = tail(fields[[1]],n=1)
  popsymbols[i] = pop
  if (pop == 'N') {p = "#698B22"}
  else if (pop == 'A1') {p = "#CD2626"}
  else if (pop == 'A2') {p = "#EEB422"}
  else if (pop == 'U') {p = "#1874CD"}
  colsymbols[i] = p
}
#head(cbind(rownames(mydata),popsymbols,colsymbols))

# Getting number of SNPs per chromosome

chr1 <- df %>% filter(chrom == "OphioH327chr_1")
rep1 <- length(chr1$num)
chr2 <- df %>% filter(chrom == "OphioH327chr_2")
rep2 <- length(chr2$num)
chr3 <- df %>% filter(chrom == "OphioH327chr_3")
rep3 <- length(chr3$num)
chr4 <- df %>% filter(chrom == "OphioH327chr_4")
rep4 <- length(chr4$num)
chr5 <- df %>% filter(chrom == "OphioH327chr_5")
rep5 <- length(chr5$num)
chr6 <- df %>% filter(chrom == "OphioH327chr_6")
rep6 <- length(chr6$num)
chr7 <- df %>% filter(chrom == "OphioH327chr_7")
rep7 <- length(chr7$num)
chr8 <- df %>% filter(chrom == "OphioH327chr_8")
rep8 <- length(chr8$num)

# Chromosome colors

chromosomes <- c(rep("darkgrey",rep1),rep("lightgrey",rep2),rep("darkgrey",rep3),rep("lightgrey",rep4),rep("darkgrey",rep5),rep("lightgrey",rep6),rep("darkgrey",rep7),rep("lightgrey",rep8))

# Setting palette color

my_palette <- c("black","#E69F00")


# Clustering and plotting genotypes

pdf("plotGenotypesHeatmap.pdf",w = 12, h = 14)

t <- heatmap.2(mydata,
          density.info = "none",
          trace = "none",
          dendrogram = "row",
          distfun = dist,
          hclustfun = hclust,
          Colv="NA",
          key=FALSE,
          col = my_palette,
          RowSideColors=colsymbols,
          ColSideColors=chromosomes,
          labCol = FALSE)

legend("topright",      # location of the legend on the heatmap plot
       legend = c("Major", "Minor"), # category labels
       fill = c("black", "#E69F00"),  # color key
       border = "white",
       bty = "n",
       horiz = TRUE)

dev.off()


# Printing order of strains from hclust
mydata_ordered <- mydata[rev(t$rowInd), t$colInd]
row_order <- rownames(mydata_ordered)
ro <- as.data.frame(row_order)
write.csv(ro, "strains_hclust_order.txt", row.names = FALSE, quote = FALSE)


