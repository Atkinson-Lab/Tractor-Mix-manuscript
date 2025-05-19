#!/bin/bash
#SBATCH --job-name=topmed-genotype
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --time=20:00:00
#SBATCH -p short
#SBATCH --mail-user=alejandra.vergara@ndph.ox.ac.uk


# Directories
base_dir=/well/emberson/projects/mcps/data/genetics_regeneron/freeze_150k/data/imputation
info_dir=/well/emberson/users/pom143/projects/LAI_TractorMix
input_dir=$info_dir/input-files
output_dir=$info_dir/geno-files/mcps
plink2=/well/emberson/shared/software/plink2/plink2



# Create a text file to keep variants with INFO score > 0.8
#awk '$16 > 0.8 {print "chr"$1}' $base_dir/MCPS_Freeze_150.GT_hg38.pVCF.rgcpid.QC2.TOPMED_dosages.snpstats.txt > $chr_out_auto/info8variants.txt

phenotype_file="$input_dir/phenotypes-ancestry-2.5P-ND-NP.txt"

    # Ensure phenotype files are created correctly
awk 'BEGIN {print "FID", "IID", "BMI", "AGE", "SEX", "Score_AMR"} {print $2, $2, $3, $4, $5, $25}' "$input_dir/extracted-phenotypes-all-2.5perc-bmi_ncj.txt" > "$phenotype_file"
    


# Run PLINK2 for each suffix
$plink2 --pgen $base_dir/MCPS_Freeze_150.GT_hg38.pVCF.rgcpid.QC2.TOPMED_dosages.HDS.pgen \
        --pvar $base_dir/MCPS_Freeze_150.GT_hg38.pVCF.rgcpid.QC2.TOPMED_dosages.HDS.pvar \
        --psam $base_dir/MCPS_Freeze_150.GT_hg38.pVCF.rgcpid.QC2.TOPMED_dosages.HDS.COLLAB.psam \
        --extract $output_dir/info8variants.txt \
        --keep "$phenotype_file" \
        --maf 0.005 \
        --hwe 1e-6 \
        --make-bed \
        --out "$output_dir/topmed-imputed-genotypes-2.5P-ND-NP"


awk 'BEGIN {print "FID", "IID", "AGE", "SEX", "Score_AMR"} {print $1, $2, $4, $5, $6}' "$phenotype_file" > "$input_dir/covariates_2.5P-ND-NP.cov"

awk 'BEGIN {print "FID", "IID", "BMI"} {print $1, $2, $3}' "$phenotype_file" > "$input_dir/phenotype_2.5P-ND-NP.pheno"


cov_file="$input_dir/covariates_2.5P-ND-NP.cov"
pheno_file="$input_dir/phenotype_2.5P-ND-NP.pheno"

# Convert lowercase suffix to uppercase P for bed, bim, and fam files
bed_suffix=$(echo "$suffix" | sed 's/p/P-ND/')

bed_file="$output_dir/topmed-imputed-genotypes-2.5P-ND-NP.bed"
bim_file="$output_dir/topmed-imputed-genotypes-2.5P-ND-NP.bim"
fam_file="$output_dir/topmed-imputed-genotypes-2.5P-ND-NP.fam"



# Run PLINK2
$plink2 --bed "$bed_file" \
    --bim "$bim_file" \
    --fam "$fam_file" \
    --covar "$cov_file" \
    --pheno "$pheno_file" \
    --glm \
    --threads 16 \
    --out "$output_dir/plink-2.5p-ND-NP"
