# Example bash commands:

### Sequencing & quality control
```
Fastqc *.fastq
```
### Read cleaning
```
parallel -j 8 java -jar -Xmx4g /prg/trimmomatic/0.33/trimmomatic-0.33.jar PE -phred33 {}_R1.fastq.gz {}_R2.fastq.gz {}_outR1P.fastq {}_outR1U.fastq {}_outR2P.fastq {}_outR2U.fastq ILLUMINACLIP:adaptest.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:3:24 MINLEN:36 ::: $(ls -1 *_R1.fastq.gz | sed 's/_R1.fastq.gz//')
```
### Read mapping
```
bwa mem Ophnu1_AssemblyScaffolds_Repeatmasked.fasta HI.4319.008.Index_15.US137_outR1P.fastq HI.4319.008.Index_15.US137_outR2P.fastq > HI.4257.001.BioOHT_97.US137_aln-pe.sam
```
### SNP and genotype calling
##### Removing duplicate reads
```
samtools view -Sb -o HI.4257.001.BioOHT_97.US137_aln-pe.bam HI.4257.001.BioOHT_97.US137_aln-pe.sam
samtools sort -o HI.4257.001.BioOHT_97.US137_aln-pe.sorted.bam HI.4257.001.BioOHT_97.US137_aln-pe.bam
samtools index HI.4257.001.BioOHT_97.US137_aln-pe.sorted.bam
java -jar picard.jar MarkDuplicates I=HI.4257.001.BioOHT_97.US137_aln-pe.sorted.bam O=HI.4257.001.BioOHT_97.US137_rmdup.bam M=HI.4257.001.BioOHT_97.US137_picard.metrics REMOVE_DUPLICATES=true
```
##### Realigning reads
```
samtools calmd -ubAr -C50 HI.4257.001.BioOHT_97.US137_rmdup.bam Ophnu1_AssemblyScaffolds_Repeatmasked.fasta > HI.4257.001.BioOHT_97.US137_calmd.bam
samtools index HI.4257.001.BioOHT_97.US137_calmd.bam
```
##### Calling genotypes (samtools v1.3; bcftools v1.4.1)
```
samtools mpileup -b bams_list_calmd -C50 -d 100000 -E -f Ophnu1_AssemblyScaffolds_Repeatmasked.fasta -q 4 -g -t DP,SP,AD,ADF,ADR,INFO/AD | bcftools call -m -f gq --ploidy 1 -o SNP_samtools_allsites.vcf -O v
```
##### Filtering SNPs
```
less SNP_samtools_allsites.vcf | grep -v "^#" | awk '{print $1"\t"$2"\t"$2"\tsnp_"NR}' > SNP_samtools_allsites_annotations
bgzip SNP_samtools_allsites_annotations
tabix -s 1 -b 2 -e 3 SNP_samtools_allsites_annotations.gz

bcftools view SNP_samtools_allsites.vcf | vcf-annotate -f +/d=10 --fill-HWE --fill-type -n -a SNP_samtools_allsites_annotations.gz -c CHROM,FROM,TO,INFO/SNP_ID -d key=INFO,ID=SNP_ID,Number=1,Type=Integer,Description='SnpList' > SNP_samtools_allsites_annotated.vcf

vcftools --vcf SNP_samtools_allsites_annotated.vcf --remove-indels --remove-filtered-all --minGQ 20 --minDP 2 --mac 1 --recode --recode-INFO-all --out SNP_samtools_SNP_filter

vcftools --vcf SNP_samtools_SNP_filter.recode.vcf --exclude-bed OphioH327_noSSR.bed --recode --recode-INFO-all --out SNP_samtools_SNP_filtered_noTE

vcftools --vcf SNP_samtools_SNP_filtered_noTE.recode.vcf --max-missing 0.5 --recode --recode-INFO-all --out SNP_samtools_SNP_filtered_noTE_miss

vcftools --vcf SNP_samtools_SNP_filtered_noTE_miss.recode.vcf --thin 273 --recode --recode-INFO-all --out SNP_samtools_SNP_filtered_noTE_miss_100000
```
##### Filtering vcf with all positions
```
bcftools view SNP_samtools_allsites.vcf | vcf-annotate -f +/d=10/-a --fill-HWE --fill-type -n -a SNP_samtools_allsites_annotations.gz -c CHROM,FROM,TO,INFO/SNP_ID -d key=INFO,ID=SNP_ID,Number=1,Type=Integer,Description='SnpList' > SNP_samtools_allsites_annotated_all.vcf

vcftools --vcf SNP_samtools_allsites_annotated_all.vcf --remove-indels --minGQ 20 --minDP 2 --remove-filtered-geno-all --recode --recode-INFO-all --out SNP_samtools_genotypes_filtered
   
vcftools --vcf SNP_samtools_genotypes_filtered.recode.vcf --exclude-bed OphioH327_noSSR.bed --max-missing 0.5 --recode --recode-INFO-all --out SNP_samtools_genotypes_filtered_noTE_miss
```
### Whole-genome alignments with *O. himal-ulmi*
##### Trimming reads with Trimmomatic
```
java -jar -Xmx4g trimmomatic-0.33.jar PE -phred33 HI.4319.008.Index_2.HP32_R1.fastq.gz HI.4319.008.Index_2.HP32_R2.fastq.gz ./seqtrimmed/HI.4319.008.Index_2.HP32_outR1P.fastq ./seqtrimmed/HI.4319.008.Index_2.HP32_outR1U.fastq ./seqtrimmed/HI.4319.008.Index_2.HP32_outR2P.fastq ./seqtrimmed/HI.4319.008.Index_2.HP32_outR2U.fastq ILLUMINACLIP:TruSeq3-PE-ophios-pauline.fa:6:20:10 MINLEN:21
```
##### Read merging with bbmerge
```
bbmerge.sh -t=8 -in1=HI.4319.008.Index_2.HP32_outR1P.fastq -in2=HI.4319.008.Index_2.HP32_outR1P.fastq -out=HI.4319.008.Index_2.HP32_merged.fastq -outu1=HI.4319.008.Index_2.HP32_R1um.fastq -outu2=HI.4319.008.Index_2.HP32_R2um.fastq
```
##### Assembly with SPAdes
```
spades.py -k 21,33,55,77,99 --careful --pe1-1 HI.4319.008.Index_2.HP32_R1um.fastq --pe1-2 HI.4319.008.Index_2.HP32_R2um.fastq --s1 HI.4319.008.Index_2.HP32_merged.fastq -o spa_HP32 -t 10 -m 200
```
##### Genome alignment and SNP calling with Mauve
```
progressiveMauve --seed-family --weight=2000 --output=himalH327.xmfa --backbone-output=himalH327.backbone Ophnu1_AssemblyScaffolds_Repeatmasked.fasta HP30_bwa_calmd.fasta HP31_bwa_calmd.fasta HP32_bwa_calmd.fasta
java -cp Mauve.jar org.gel.mauve.analysis.SnpExporter -f himalH327.xmfa -o himalH327.snps
```
### Mitochondrial DNA
##### Changing variants to fasta alignment excluding positions with any missing data
```
python writeVcfTofasta.py
```
##### Maximum likelihood phylogeny
```
raxmlHPC-SSE3 -s SNP_samtools_mtDNA.fasta -n mtDNA -m GTRCAT -V -p 12345 -f a -# 100 -x 123459 -o HP32_H,HP30_H,HP31_H
```
##### Collapsing branches with <50 bootstrap support
```
java -jar TreeCollapseCL4.jar -f RAxML_bipartitions.mtDNA -b 50
```
### Genetic diversity and divergence
##### Calculating nucleotide diversity and Tajima’s D with ANGSD
```
angsd -bam bamlist_ame -gl 1 -ref Ophnu1_AssemblyScaffolds_Repeatmasked_HM.fasta -anc himalH327_HP30_masked.fasta -dosaf 1 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -baq 1 -C 50 -minMapQ 1 -P 4 -out ame
misc/realSFS ame.saf.idx -maxIter 100 -P 4 >small_ame.sfs
angsd -bam bamlist_ame -out ame -doThetas 1 -doSaf 1 -pest small_ame.sfs -anc himalH327_HP30_masked.fasta -ref Ophnu1_AssemblyScaffolds_Repeatmasked_HM.fasta -GL 1 -P 4 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -baq 1 -C 50 -minMapQ 1
```
*Getting per chromosome estimate*
```
misc/thetaStat do_stat ame.thetas.idx
```
*Getting per site estimates*
```
misc/thetaStat do_stat ame.thetas.idx -win 1 -step 1 -outnames ame.thetasSite.gz
```
*Getting per window estimates*
```
misc/thetaStat do_stat ame.thetas.idx -win 50000 -step 50000 -outnames ame.thetasWindow.gz
```
##### Calculating per-site Dxy with ANGSD
*Running angsd to obtain variable positions in all samples*
```
angsd -bam bamlist -doMaf 3 -ref Ophnu1_AssemblyScaffolds_Repeatmasked_HM.fasta -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -baq 1 -C 50 -minMapQ 1 -P 4 -out tot -doMajorMinor 1 -GL 1 -skipTriallelic 1 -SNP_pval 1.0e-3
```
*Selecting sites of interest*
```
less tot.mafs.gz | tail -n+2 | cut -f 1-2 > sites.txt
angsd sites index sites.txt
```
*Running angsd in populations to find maf (minor allele frequency)*
```
angsd -bam bamlist_ame -doMaf 3 -ref Ophnu1_AssemblyScaffolds_Repeatmasked_HM.fasta -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -baq 1 -C 50 -minMapQ 1 -P 4 -out ame 	-doMajorMinor 5 -anc himalH327_HP30_masked.fasta -GL 1 -sites sites.txt
gunzip -c ame.mafs.gz > ame.mafs
```
*Getting persite Dxy*
```
Rscript ngsTools/ngsPopGen/scripts/calcDxy.R -p ame1.mafs -q ame2.mafs -t 29426211 > ame1_ame2.dxy.global
less Dxy_persite.txt | awk '{if ($3!=0) print $0}' > ame1_ame2.Dxy_persite.txt
```
##### Calculating Fst with ANGSD
```
misc/realSFS -P 4 ame1.saf.idx ame2.saf.idx > ame1.ame2.ml
misc/realSFS fst index ame1.saf.idx ame2.saf.idx -sfs ame1.ame2.ml -fstout ame1.ame2
```
*Getting global estimate*
```
misc/realSFS fst stats ame1.ame2.fst.idx 2> ame1.ame2.fst.global
```
*Getting per window estimate*
```
misc/realSFS fst stats2 ame1.ame2.fst.idx -win 50000 -step 50000 -type 2 > ame1.ame2.fst.window
```
*Getting per site estimate*
```
misc/realSFS fst stats2 ame1.ame2.fst.idx -win 1 -step 1 -type 2 > ame1.ame2.site.fst.window
```
##### Converting genotypes to major and minor alleles
*Filtering SNPs with missing information*
```
vcftools --gzvcf SNP_samtools_SNP_filtered_noTE_miss.recode.vcf.gz --max-missing 1 --recode --recode-INFO-all --out SNP_samtools_SNP_filtered_noTE_miss1
```
*Annotating SNPs with AF field*
```
bcftools +fill-tags SNP_samtools_SNP_filtered_noTE_miss1.recode.vcf -Ov -o SNP_samtools_SNP_filtered_noTE_miss1.recode.AF.vcf -- -t AF
```
*Extracting table with AF and genotypes*
```
bcftools view -m2 -M2 -v snps SNP_samtools_SNP_filtered_noTE_miss1.recode.AF.vcf -Ov | bcftools query -f '%CHROM\t%POS\t%INFO/AF[\t%GT]\n' > genotypes1.tab
```
*Converting table with genotypes to major and minor alleles*
```
python convert2MajorMinor.py
```
### Phylogenetic discordances across the genome
##### Discordant topologies across the genome with Twisst
```
python twisst.py -t raxml_trees.trees -w output.weights.csv.gz -g A1 -g A2 -g N -g U -g H --method complete --groupsFile grouping.txt --outputTopos topologies.trees
```
### Functional analysis of introgressed regions
##### Blasting H327 proteins against nr database with diamond
```
diamond makedb --in ./DATA/nr/nr --db ./DATA/nr/nr
diamond blastp -p 18 -d ./DATA/nr/nr -q Ophnu1_GeneCatalog_proteins_20170425.aa.fasta -o Ophnu1_GeneCatalog_proteins_20170425.aa.fasta.xml -f xml --sensitive -l 1 -k 10 --salltitles
```
### Selection test
##### Estimating synonymous and nonynonymous sites in genes with mstatspop
```
mstatspop -f fasta -i Gene_1981982.fasta -o 1 -p 1 -u 1 -s 123456 -G 0 -N 2 73 21 -n fai/Gene_1981982.txt -g gff/Gene_1981982.gff synonymous Nuclear_Universal > output/Gene_1981982_syn.out'
mstatspop -f fasta -i Gene_1981982.fasta -o 1 -p 1 -u 1 -s 123456 -G 0 -N 2 73 21 -n fai/Gene_1981982.txt -g gff/Gene_1981982.gff nonsynonymous Nuclear_Universal > output/Gene_1981982_nsyn.out'
```
### Inversion between AME1 (and ULM) and AME2 (and NOV)
##### Aligning assemblies with nucmer
```
nucmer --maxgap=1000 --mincluster=100 --prefix=pair_DDS241_04_51 spades_DDS241_AME1.fasta spades_04_51_AME2.fasta
delta-filter -q -r pair_DDS241_04_51.delta > pair_DDS241_04_51_qr.filter
show-coords -rclT pair_DDS241_04_51_qr.filter > pair_DDS241_04_51_qr.coords
```

