## Tractor-Mix Simulation

We conducted *in silico* simulations to evaluate the false positive rates, statistical power, and bias of estimators in the Tractor-Mix model. The code is written in R but was executed on a high-performance computer. Therefore, the R code is named `xxx.R`, while the bash script to execute the R code is named `exec.xxx.sh`. Special attention should be given to the R code.

This README serves as a guideline to replicate the simulation results in Tractor-Mix. Please contact Taotao Tan (Taotao.Tan@bcm.edu) if anything is unclear.

&nbsp;  
&nbsp;  

### Genotype Simulation

**`Utils.R`:**  
`GenerateChild()` is used for simulating the genotype for admixture with the gene-dropping method.  
`GenerateMafs()` is used for generating diverged allele frequencies according to the Balding–Nichols model.  
`GenerateHomo()` is used for generating a homogeneous population. 
`GenerateAdm()` is used for generating an independent admixed population.  
`WriteGenotype()` is used for producing `dosage` or `hapcount` files.  
`CalcLambdaGC()` is used for calculating lambda GC, allowing for calculation with different degrees of freedom.  
`PrepareGDS()` is used for converting simulated genotype files into Plink format and computing PC-Air and PC-Relate.

**`Pedigree3_MAF.R`:**  
This code simulates allele frequencies according to the Balding–Nichols model. 

**`Pedigree3.R`:**  
I used `MakePedigree()` and `MakeGenotype()` to generate the genotype matrix for an admixed family (Pedigree is defined in Supplement 2). The outputs are `G_Tot.tsv`, `G_Afr.tsv`, `G_Eur.tsv`, and `LA.tsv`, representing total risk allele counts, African-specific allele counts, European-specific allele counts, and local ancestry matrix for AFR, respectively.  
`MakeGRM()` is used for producing the true kinship matrix according to the pedigree.  

**`Pedigree3_Rel.R`:**  
This code simulates related admixed samples using functions from `Pedigree3.R`. 

**`Pedigree3_Ind.R`:**  
This code simulates independent admixed samples.

**`create_Block_row.R`:**   
`Pedigree3_Rel.R` and `Pedigree3_Ind.R` create 50 `Blocks` (each block contains 100 samples and 1M variants). However, this data format is not computationally efficient. `create_Block_row.R` reorganizes `Blocks` into 20 `Block_rows`, with each `Block_row` containing all 5000 samples and approximately 50,000 variants. The `Block_row` are formatted as `dosage` and `hapcount`, matching the desired input format of `Tractor-Mix`. 

**`qc_Block_row.R`:**  
Some QC is performed to eliminate variants with low allele frequency (between 0.05 and 0.95).

&nbsp;  
&nbsp;  

### False Positive Simulations 

`Block 1-2` (approximately 90,000 variants) were used to compute PC, PC-Air, GRM, and PC-Relate. 

`Block 3-18` (approximately 810,000 variants) were used for conducting statistical tests. 

**`Calc_GENESIS.R`:**  
First convert data to Plink format, then use `snpgdsIBDKING` to compute the King-robust estimator. Then compute PC-Air with King-robust, and compute PC-Relate with PC-Air, using the top 2 principal components.

**`Calc_PC_GRM.R`:**  
First standardize the genotype matrix, then compute GRM according to ($GG^T / p$). `prcomp` was used for PC analysis. 

**`Fit_Null_cont.R`, `Fit_Null_dich.R`:**  
Generate continuous and binary phenotypes, then fit the 3 null models:  
* Homogeneous PC + Homogeneous GRM  
* PC-Air + PC-Relate  
* Admixture proportion + Kinship  

The intercept term of dichotomous phenotype was adjusted according to prevalence (0.2).

**`Compare_FP_cont.R`, `Compare_FP_dich.R`:**  
Run 8 different models and generate their summary statistics (see manuscript). 

**`SparseFP_cont.R`, `SparseFP_dich.R`:**   
Only test 1 model with the kinship matrix and admixture proportion. Sparse matrices were used to increase computational efficiency. 

**`Aggregate_FP_res.R`:**  
Aggregate the results from `SparseFP_cont.R` and `SparseFP_dich.R`. It computes the false positive rates under different $\alpha$ levels. 

&nbsp;  
&nbsp;  

### Statistical Power

**`Power_simu_cont.R`, `Power_simu_dich.R`:**  
These two scripts compute the *joint* p-values and effect size estimates for both the score test and the Wald test, according to the parameters in `FP_params_cont.tsv` and `FP_params_dich.tsv`. 

&nbsp;  
&nbsp;  

### Compare p-values from score and Wald tests 

**`Comp_ancP.R`:**
Generate phenotype under null and alternative hypotheses, then use the score test and Wald test to compute p-values. 

&nbsp;  
&nbsp;  

### Tractor-Mix Code

These are all the `Tractor-Mix` scripts used in simulations. These scripts are prototypes of `Tractor-Mix` and should not be used in realistic GWAS analysis. 

**`TractorMix.score.R`**  
Implements a basic score test for Tractor-Mix. It calculates a joint p-value and chi-square statistics.

**`TractorMix.ancP.score.R`**  
Implements an approach to calculate ancestry-specific p-values: ($s.e. = \sqrt{diag [Var(T)^{-1} ]}$).

**`TractorMix.sparse.score.R`**  
Implements Tractor-Mix that allows sparse GRM as inputs for the score test.

**`TractorMix.wald.R`**  
Implements the Wald test for Tractor-Mix. It directly fits the full model and produces ancestry-specific effect sizes and p-values.



