# Tractor-Mix UKBB analysis
Scripts for Tractor-Mix analysis on UKBB AFR-EUR admixture.

Inputs: 
Genotype file: `/storage/atkinson/shared_resources/datasets/UKB/QCedGenotypes/UKB-AFR-biallelicSNPs-MAF005-95miss-INFO80-HWE1e6-auto.vcf.bgz`

TGP reference panel (split by 22 chromosomes): `/storage/atkinson/shared_resources/reference/ReferencePanels/TGP/TGP_Phase3_GRCh37/FilteredVCFs/ALL.chr$CHR.phase3_shapeit2_mvncall_integrated_v5b.maf005.biallelic.unrelated.recode.vcf.gz`


&nbsp;  
&nbsp;  

### `PC_King` compute PCA and King-robust

`PC_King.sh`:  
This script performs LD prunning, then uses the pruned VCF to compute King-robust and PC.


&nbsp;  
&nbsp;  

### `ADMIXTURE` compute ADMIXTURE proportion

`integrated_call_samples_v3.20130502.ALL.panel`:  
This is the original meta-data for TGP samples - contains basic information including population and super population (continental level)

`TGP_AFR_EUR_4Plink.tsv`:  
A list of samples to be extracted from TGP. This list of samples only contains AFR and EUR.  

`TGP_processing.sh`:  
Pre-process TGP data. Merge 22 chromosomes, and extract variant list

`UKBB_processing.sh`:  
Perform LD pruning for UKBB data, and extract variant list

`intersect_vars.R`:
XXXXXX NEED TO WRITE XXXXXX

`ADMIXTURE.sh`:  


&nbsp;  
&nbsp;  


### `PHASE` phasing UKBB data  

`split_chrom.sh`:  
Split UKBB vcf file by 22 chromosomes. This helps parallelization for downstream analysis 

`run_shapeit5.sh`:  
Phasing each chromosome with Shapeit5. First use bcftools to subset TGP (reference) AFR-EUR samples, then compute allele count (AC) and indexing both TGP and UKBB data. Then perform phasing, which outputs a bcf file. Then convert bcf files to vcf.gz files. Genetics maps downloaded from [here](https://github.com/odelaneau/shapeit4/tree/master/maps)



&nbsp;  
&nbsp;  

### `LAI` local ancestry inference

`run_RFMIX.sh`:  
Perform local ancestry inference with RFmix2. Genetics maps downloaded from [here](https://github.com/odelaneau/shapeit4/tree/master/maps), and modified to match the format of RFMix2. 

`TGP_AFR_EUR.tsv`:  
The super population level label for TGP reference panel. 








