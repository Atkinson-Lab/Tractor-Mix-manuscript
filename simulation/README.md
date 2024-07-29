## Tractor-Mix Simulation

We conducted *in silico* simulations to evaluate the false positive rates, statistical power, bias of estimators of the Tractor-Mix model. The code is written in R, but was executed with high performance computer. Therefore, the R code is named as `xxx.R`, while the bash script to execute R code is named as `exec.xxx.sh`. One should pay more attention to the R code.

This README served as a guideline to replicate the simulation results in Tractor-Mix. Please contact Taotao Tan (Taotao.Tan@bcm.edu) if anything is unclear

&nbsp;  
&nbsp;  

### Genotype Simulation

**`Utils.R`:**  
`GenerateChild()` is used for simulating the genotype for admixture with gene dropping method.  
`GenerateMafs()` is used for generating diverged allele frequencies according to Balding–Nichols model.  
`GenerateHomo()` is used for generating homogeneous population. 
`GenerateAdm()` is used for generating independent admixed population.  
`WriteGenotype()` is used for producing `dosage` or `hapcount` file.  
`CalcLambdaGC()` is used for calculating lambda GC. It allows for calculating lambda GC with different degree of freedom.  
`PrepareGDS()`  is used for converting simulated genotype files into Plink format, and compute PC-Air and PC-Relate 

**`Pedigree3_MAF.R`:**  
This is the code to simulate allele frequencies according to Balding–Nichols model. 

**`Pedigree3.R`:**  
I used `MakePedigree()` and `MakeGenotype()` for generating the genotype matrix for an admixed family (Pedigree is defined in Supp 2). The outputs are `G_Tot.tsv`, `G_Afr.tsv`, `G_Eur.tsv`, `LA.tsv`, stands for total risk allele counts, Afr specific allele counts, Eur specific allele counts, local ancestry matrix for AFR.  
`MakeGRM()` is used for producing the true kinship matrix according to the pedigree.  

**`Pedigree3_Rel.R`:**  
This is the code to simulate related admixed samples. It used functions from `Pedigree3.R`. 

**`Pedigree3_Ind.R`:**  
This is the code to simulate independent admixed samples.

**`create_Block_row.R`:**   
`Pedigree3_Rel.R` and `Pedigree3_Ind.R` creates 50 `Block`s (each block contains 100 samples and 1M variants). However, this data format is not computationally efficient. The `create_Block_row.R` re-organize `Block`s into 20 `Block_row`s, with. each `Block_row` contains all 5000 samples and approximately 50,000 variants. The `Block_row` are formatted as `dosage` and `hapcount`, matching the desired input format of `Tractor-Mix`. 

**`qc_Block_row.R`:**  
Some QC is performed to eliminate variants with low AF (between 0.05 to 0.95).

&nbsp;  
&nbsp;  

### False Positive simulations 

`Block 1-2` (approximately 90,000 variants) were used to compute PC, PC-Air, GRM, PC-Relate 

`Block 3-18` (approximately 810,000 variants)  were used for conducting statistical tests. 

**`Calc_GENESIS.R`:**  
First convert data to Plink format, then use `snpgdsIBDKING` to compute the King-robust. Then compute PC-Air with King-robust, and compute PC-Relate with PC-Air, using top 2 principal components.

**`Calc_PC_GRM.R`:**  
First standardize the genotype matrix, then compute GRM according to $GG^T / p$. `prcomp` was used for PC analysis. 

**`Fit_Null_cont.R`, `Fit_Null_dich.R`:**  
Generate continuous and binary phenotypes. Then fit the 3 null models:  
* Homogeneous PC + Homogeneous GRM;  
* PC-Air + PC-Relate;  
* Admixture proportion + Kinship;  
The intercept term of dichotomous phenotype were adjusted according to prevalence (0.2)

**`Compare_FP_cont.R`, `Compare_FP_dich.R`:**  
Run 8 different models, and generate their summary statistics (See manuscript). 

**`SparseFP_cont.R`, `SparseFP_dich.R`:**  
Only test 1 model with the kinship matrix and admixture proportion. I used sparse matrix to increase computational efficiency. 


&nbsp;  
&nbsp;  

### Statistical Power

**`Power_simu_cont.R`, `Power_simu_dich.R`:**  
Those two scripts computed the p-value and effect size estimates for both score test and Wald test, according to the parameters in `Cont_FP_res.tsv, Dich_FP_res.tsv`. 



