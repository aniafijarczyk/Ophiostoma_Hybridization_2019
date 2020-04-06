#### Custom code used in the paper: [Hybridization and introgression drives genome evolution of the Dutch Elm Disease pathogens](https://www.nature.com/articles/s41559-020-1133-6?proof=trueNov) 


Open source software
--------------------

All open source software used in the study is listed in _Open_source_software.md_


Bash commands
-------------

Bash commands used to run the analyses are in the file _Bash_commands.md_


Custom code
-----------

All custom scripts are organized in four folders according to the appearance in the Results section


### 01_Population_structure_and_cryptic_lineages:
**DAPC_analyses:**
- plotDAPC_SNP.R  # plotting DAPC analysis based on SNP variants
- plotDAPC_CNV.R  # plotting DAPC analysis based on CNV variants

**Figure1_plots:**
- simulate_point_locations.R	# simulating coordinates of point labels
- map_geo.R  			# plotting a map with strains
- plotGlobe.R			# plotting picture of a globe
- plotTree.py  			# plotting NJ phylogeny
- plotStructure.R  		# plotting structure analysis
- PCA_Plot.R  			# plotting PCA results

**Mitochondrial_DNA**
- writeVcfToFasta.py	# conversion of vcf to fasta
- drawTree.py		# plotting mtDNA tree

**Strcuture_and_map_O_novo_ulmi**
- simulate_point_locations.R    # simulating coordinates of point labels
- plotStructureMap.R		# plotting structure results for K = 3 on a map
- plotStructure_novoulmi.R	# plotting strcuture barplots for different K

### 02_Diversity_and_genomic_signatures_of_introgression
**Diversity_and_genotype_heatmap**
- combineDxy.py			# getting Dxy per window estimates
- convert2MajorMinor.py		# converting genotypes into major-minor allele
- plotGenotypesHeatMap.R	# plottig heatmap of genotypes
- plotDiversityDist.R		# plotting distribution of diversity in 50 kb windows
- calculateCorrelations.R	# plotting diversity/divergence correlations

**Diversity_plots_in_windows**
- plotAngsdWindows_short_Pi_Dxy.R	# plotting tree weighting, Pi and Dxy, in windows across the genome
- plotAngsdWindows_short_Fst_Taj.R	# plotting tree weighting, Tajima's D and Fst in windows across the genome

**Figure2_plots**
- plotMatrix.R		# plotting a matrix of global and per-window distances and diversity in Ophiostoma lineages
- plot_Dstatistic.R	# dotplot with standard errors of D-statistic
- mk_circos.sh		# script for plotting circos chronogram

**Figure3_plots**
- plotTargetsAndDiversity.R	# plotting ULM ancestry, diversity within, divergence between lineages and targets of positive selection across the genome

**Gene_flow_between_Ophiostoma_lineages**
- runUnbiasedD_Outgr.py		# estimating D-statistic

**High_diversity_regions**
- bootstrapPi.py	# detecting high diversity regions across the genome

**MAT_loci**
- plotDepth_heatmap.R	# plotting coverage and admixture around mating types
- plotTree_MAT.py	# plotting phylogenetic trees around mating types

**PCA_in_windows**
- runPCA.R	# running PCA analysis for a window and extracting distance of each strain with ULM

**Structure_in_windows_plots**
- mk_circos.sh	# script for plotting circos chronogram with structure in windows results

### 03_Genetic_basis_of_lineage_divergence
**Identification_of_ULM_ancestry**
- filter_ULMLike.py	# Identifying strains with ULM ancestry from tree topology
- filterFasta.py	# selecting strains from fasta alignments

**Selection_test**
- SnIPRE_running_script_Ophiostoma.R	# running snipre program to detect genes under selection
- plotSnipreResults.R			# plotting snipre results

**Genetic diversity and divergence**
- diversity_divergence_stats_popgenome.R	# calculating diversity/divergence in genes

### 04_Phenotypic_differences_between_lineages
**Processing pictures and extracting growth information**
- ophiostoma_phenotyping_picture_analysis_function.R 	# a function to process images
- ophiostoma_phenotyping_launch.R 			# launching a function for each picture

**Modelling and plotting isolate growth in different conditions and virulence on apples**
- growth_virulence_modeling_Fig4_S17_S18_script.R	# plotting growth & virulence results


