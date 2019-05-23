'''
This scipt reads vcf with SNPs, genome in fasta file with ancestral state,
and a file with population info, and calculates unbiased D (Patterson's D) 
for different combinations of 4 taxa

USAGE:

python runUnbiasedD_Outgr.py

requirements:       sample.recode.vcf
                    Species_distinction_PCA.txt
                    himalH327_HP30_OphioH327chr_[1-8].fasta
'''

import numpy as np
import pandas as pd
import allel; print('scikit-allel', allel.__version__)
from Bio import AlignIO


print("Reading vcf, only biallelic SNPs")

vcf = "sample.recode.vcf"
calltest = allel.read_vcf(vcf,alt_number=1,numbers={'calldata/GT': 1})
haplotypes = allel.HaplotypeArray(calltest['calldata/GT'])
chromosomes = calltest['variants/CHROM']
refs = calltest['variants/REF']
alts = calltest['variants/ALT']
variants = calltest['variants/POS']
Nvariants = len(variants)
vars = range(Nvariants)


print("Setting populations")

samples_fn = "Species_distinction_PCA.txt"
samplesPops = pd.read_csv(samples_fn,sep='\t',header=0,index_col = 'SpeciesID')
samplesDict = samplesPops.to_dict('index')
samples = calltest['samples']
samples_subset = pd.DataFrame({'FileName':samples[:]})
groups = [samplesDict[ele]['GroupPCA1'] for ele in samples_subset['FileName']]
samples_subset['Group'] = groups

subpops = {
    'all': list(range(len(samples_subset))),
    'nov': samples_subset[samples_subset.Group == 'N'].index.tolist(),
    'ame1': samples_subset[samples_subset.Group == 'A1'].index.tolist(),
    'ame2': samples_subset[samples_subset.Group == 'A2'].index.tolist(),
    'ulm': samples_subset[samples_subset.Group == 'U'].index.tolist(),
    'him': samples_subset[samples_subset.Group == 'H'].index.tolist(),
    'ame1H': samples_subset[samples_subset.Group == 'A1Hyb1'].index.tolist(),
    'ame2H': samples_subset[samples_subset.Group == 'A2Hyb1'].index.tolist(),
    'novH': samples_subset[samples_subset.Group == 'NHyb1'].index.tolist(),
    'ulmH': samples_subset[samples_subset.Group == 'UHyb1'].index.tolist(),    
    'ame': samples_subset[samples_subset.Group == 'A1'].index.tolist() + samples_subset[samples_subset.Group == 'A2'].index.tolist()
}


print("Reading outgroup")

CH = {}
for chromo in range(1,9):
    print(chromo)
    try:
        myFile = AlignIO.read("himalH327_HP30_OphioH327chr_"+str(chromo)+".fasta","fasta")
    except IOError:
        continue
    else:
        nsites = len(myFile[0])
        chrID = "OphioH327chr_"+str(chromo)
        P = {}
        for i in range(nsites):
            pos = i+1
            alleles = list(myFile[:,i])
            CH[chrID+"."+str(pos)] = [x.upper() for x in alleles]


print("Getting allele frequencies for outgroup")

HimListFreq = []
HimListCount = []
for ele in range(Nvariants):
    ind = chromosomes[ele]+"."+str(variants[ele])
    alleles = list(refs[ele])+list(alts[ele])
    alleles2 = [i for i in alleles if i]
    himalAlleles = [x for x in CH[ind] if x != 'N']
    # the ancestral variants cannot be missing and the total number of segregating variants with ancestral one is not more than 2
    test = len(set(list(alleles2) + himalAlleles))
    #b = [x for x in himalAlleles if x not in alleles]
    if himalAlleles and (test <= 2):
        refHimalFreq = himalAlleles.count(refs[ele])/float(len(himalAlleles))
        refHimalCount = np.array([himalAlleles.count(refs[ele]),himalAlleles.count(alts[ele][0])],dtype=np.int32)
    else:
        refHimalFreq = np.nan
        refHimalCount = np.array([np.nan, np.nan])
    HimListFreq.append(refHimalFreq)
    HimListCount.append(refHimalCount)
HimArray = np.array(HimListFreq)
HimCArray = np.array(HimListCount,dtype=np.int32)


print("Choosing tests")

Tests = ['out_ulm_ame_nov',
         'out_ulm_ame1_ame2',
         'out_ulm_ame1_nov',
         'out_ulm_ame2_nov',
         'out_nov_ame1_ame2',
         'ulm_nov_ame1_ame2',
         'out_ame2_ame1H_ame1',
         'out_nov_ame1H_ame1',
         'out_ulm_ame1H_ame1',
         'out_ame1_ame2H_ame2',
         'out_nov_ame2H_ame2',
         'out_ulm_ame2H_ame2',
         'out_ame1_novH_nov',
         'out_ame2_novH_nov',
         'out_ulm_novH_nov',  
         'out_ame1_ulmH_ulm',
         'out_ame2_ulmH_ulm',
         'out_nov_ulmH_ulm']



wh = open('runUnbiasedD_Outgr_parents_genome.out','w')
wh.write('Test\tNsites\td\tse\tZ\tNblocks\tNsitesInBlock\n')

for test in Tests:
    
    popW = test.split('_')[0]
    popX = test.split('_')[1]
    popY = test.split('_')[2]
    popZ = test.split('_')[3]
    testID = '_'.join([popW,popX,popY,popZ])
    blen = 200
    
    print("Getting allele frequencies for populations")
    
    pops = {'ame1':0,'ame2':0,'ame':0,'nov':0,'ulm':0,'ame1H':0, 'ame2H':0, 'novH':0, 'ulmH':0,}
    New_array = {}
    for pop in pops.keys():
        haploid_array = haplotypes.subset(vars,subpops[pop])
        allel_count = haploid_array.count_alleles(max_allele=1)
        New_array[pop] = allel_count
        allel_freq = allel_count.to_frequencies()[:, 0]
        pops[pop] = allel_freq
    pops['out'] = HimArray
    New_array['out'] = HimCArray

    print("Filtering variants which are fixed in a pair of populations")

    nvlist = []
    for i in vars:
        b=0
        if (pops[popY][i] + pops[popZ][i] == 0) | (pops[popY][i] + pops[popZ][i] == 2):
            b+=0
        else:
            b+=1
        if (pops[popX][i] + pops[popW][i] == 0) |  (pops[popX][i] + pops[popW][i] == 2):
            b+=0
        else:
            b+=1

        if b==2:
            nvlist.append(i)

    print("Calculating D")
    
    w = New_array[popW][nvlist]
    x = New_array[popX][nvlist]
    y = New_array[popY][nvlist]
    z = New_array[popZ][nvlist]
    ud = allel.average_patterson_d(w, x, y, z, blen)
    f = "D("+popW+","+popX+","+popY+","+popZ+")", str(len(nvlist)), str(ud[0]), str(ud[1]), str(ud[2]), str(len(ud[3])), str(blen)
    t = w,x,y,z
    wh.write('\t'.join(f)+'\n')


    print("Test\tD("+popW+","+popX+","+popY+","+popZ+")")
    print("N_sites\t"+str(len(nvlist)))
    print("N_blocks\t"+str(len(ud[3])))
    print("N_variants_in_block\t"+str(blen))
    print("d\t"+str(ud[0]))
    print("se\t"+str(ud[1]))
    print("Z\t"+str(ud[2]))
    print("\n")
    #print(ud)

wh.flush()
wh.close()



