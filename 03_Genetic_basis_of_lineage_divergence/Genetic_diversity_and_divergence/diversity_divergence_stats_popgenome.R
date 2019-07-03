library(PopGenome)

GENOME.class=readData("FASTA_examples")
GENOME.class

populations=read.table("populations_diversity.txt")

novpop=subset(populations, V2=="_N")
novpop=droplevels(novpop)
nov=as.vector(novpop$V1)

ulmpop=subset(populations, V2=="_U")
ulmpop=droplevels(ulmpop)
ulm=as.vector(ulmpop$V1)

ame1pop=subset(populations, V2=="A1")
ame1pop=droplevels(ame1pop)
ame1=as.vector(ame1pop$V1)

ame2pop=subset(populations, V2=="A2")
ame2pop=droplevels(ame2pop)
ame2=as.vector(ame2pop$V1)

GENOME.class <- set.populations(GENOME.class,list(nov,ame1, ame2, ulm), diploid = F)

### Neutrality stats
GENOME.class <- neutrality.stats(GENOME.class, detail = T)
neut=get.neutrality(GENOME.class)[[2]]
get.neutrality(GENOME.class)[[4]]
GENOME.class@Tajima.D

### MK test
GENOME.class <- MKT(GENOME.class, do.fisher.test = T)
mkt_data=as.data.frame(get.MKT(GENOME.class), row.names = NULL)

### Diversities
GENOME.class <- diversity.stats.between(GENOME.class)
GENOME.class <- diversity.stats(GENOME.class)

