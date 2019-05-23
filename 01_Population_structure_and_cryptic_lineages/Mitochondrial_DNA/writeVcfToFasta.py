'''
This script read vcf with SNPs, and a file with population info,
converts genotypes to fasta alignment and filters positions with any missing data.

USAGE:
    
python writeVcfToFasta.py

requirements:   SNP_samtools_mtDNA_variants_filter2_miss.recode.vcf
                Species_distinction_PCA.txt
output:         SNP_samtools_mtDNA_noMD.phy
'''



import numpy as np
import pandas as pd
import allel; print('scikit-allel', allel.__version__)
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Seq import Seq

#popTrans = {'ame1': 'A1', 'ame2':'A2', 'nov':'N', 'ulm':'U','him':'H'}

himal = ['HI.4257.001.BioOHT_42.HP30_calmd.bam',
         'HI.4257.001.BioOHT_41.HP32_calmd.bam',
         'HI.4257.001.BioOHT_8.HP31_calmd.bam']

### Reading vcf file

vcf = "SNP_samtools_mtDNA_variants_filter2_miss.recode"
callset = allel.read_vcf(vcf+".vcf")
genotypes = allel.GenotypeChunkedArray(callset['calldata/GT'])
variants = callset['variants/POS']
samples = callset['samples']
alts = callset['variants/ALT']
refs = callset['variants/REF']
chroms = callset['variants/CHROM']

### Reading pop info

samples_fn = "Species_distinction_PCA.txt"
df = pd.read_csv(samples_fn,sep='\t',header=0)
ndf = df.set_index('SpeciesID.1')
d = ndf.to_dict()['GroupSamtools']


### Changing samples names in vcf

T = {}
for s in samples:
    news = s.replace('HI.4257.001.BioOHT_','').replace('_calmd.bam','').split('.')[1]
    T[s] = news


### Converting genotypes to fasta
    
#wh = open('SNP_samtools_mtDNA.fasta','w')
S = {}
for samp in range(len(samples)):
    Sek = []
    X = {}
    for snp in range(len(variants)):
        gt = genotypes[snp,samp][0]
        alt = np.append(refs[snp][0],alts[snp])

        lengths = [len(ele) for ele in alts[snp] if len(ele)>1]
        if lengths:
            continue
        else:
            if gt==-1:
                g = 'N'
            elif gt==0:
                g = refs[snp][0]
            elif gt==1:
                g = alts[snp][0]
            elif gt==2:
                if alts[snp][1]:
                    g = alts[snp][1]
                else:
                    g = 'N'
            elif gt==3:
                if alts[snp][2]:
                    g = alts[snp][2]
                else:
                    g = 'N'
            else:
                g = 'N'
            Sek.append(g)
    header = T[samples[samp]]+"_"+d[T[samples[samp]]]
    sek = ''.join(Sek)
    #wh.write('>'+header+'\n')
    #wh.write(sek+'\n')
    S[header] = Sek
#wh.flush()
#wh.close()


### Filtering positions with missing data

df = pd.DataFrame(S)
dt = df.transpose()
num_sites = len(list(dt))
dn = dt.apply(pd.value_counts)
N_sites = [] #positions where there is at least one N site
few_N_sites = [] #positions where there are less than 10% missing data

for pos in range(num_sites):
    ens = dn[pos]['N']
    ens_perc = float(ens)/float(num_sites)
    if ens > 0:
        N_sites.append(pos)
    if ens_perc > 1:
        few_N_sites.append(pos)

dt2 = dt
ds1 = dt.drop(N_sites, axis=1)

print("Number of SNPs remaining: "+str(len(list(ds1)))+"\n")

df_new = ds1.transpose()
ResList = []
speciesCols = list(df_new)
for sp in speciesCols:
    sek = ''.join(df_new[sp])
    record = SeqRecord(Seq(sek),id=sp)
    ResList.append(record)


# Writing dataframe to a file

hand = open("./SNP_samtools_mtDNA_noMD.phy", "w")
SeqIO.write(ResList, hand, "phylip")
hand.flush()
hand.close()


