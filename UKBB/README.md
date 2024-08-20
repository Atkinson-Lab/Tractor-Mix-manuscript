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
Intersect shared variants in TGP and UKBB

`ADMIXTURE.sh`:  
Extract shared variants for UKBB and TGP, then merge them and convert them to vcf files, then run ADMIXTURE


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


&nbsp;  
&nbsp;  

### `TRACTOR` Tractor-Mix run 

&nbsp;  

##### Prepare ancestry specific dosage files (for Tractor-Mix)

`run_ExtractTracts.sh`:  
Allocate alleles to each local ancestry. It takes the `msp` files and phased `vcf` files as input, and produce `anc*.dosage.txt` and `anc*.hapcount.txt`. 

&nbsp;  

##### Prepare bed files for each chromosomes (for GMMAT)

`run_split_convert_bed.sh`
Convert the original VCF to bim/bam/fam, then split the file by chromosome. The splited chromosome files are used for GMMAT run

&nbsp;  

##### Prepare PC-Air and PC-Relate  

`Make_King.sh`:  
First make a `king` file from the pruned UKBB vcf file; then convert vcf file to bed/bim/fam file 

`make_PC_GRM_from_GENESIS.R`, `run_make_PC_GRM_from_GENESIS.sh`:  
Use `bed/bim/fam` and `king` as inputs, run GENESIS to compute PC-Air and PC-Relate

&nbsp;  

##### Prepare Null models

`fit_Nulls.R`:  
Harmonize phenotype dataframe, then fit null models with PC-Air and PC-Relate for TC (total cholestral), LDL (low-density lipoprotein), SCA (sickle cell anemia)

&nbsp;  

##### Run score test 

`run_TC_LDL_SCA_GMMAT.sh`, `TC_LDL_SCA_GMMAT.R`:  
Run GMMAT for 3 phenotypes. 

`run_TC_LDL_SCA_TractorMix.sh`, `TC_LDL_SCA_TractorMix.R`:
Run Tractor-Mix for 3 phenotypes. 

`run_SCA_TractorMix_AC.sh`, `SCA_TractorMix_AC.R`: 
Run GWAS for SCA, but also count minor alleles for each variant


&nbsp;  
&nbsp;  


### Tractor-Mix source code 

`TractorMix.score.R`: flagship code 

`TractorMix.score.AC.R`: compute allele counts while conducting Score test 

`TractorMix.wald.R`: the wald test




