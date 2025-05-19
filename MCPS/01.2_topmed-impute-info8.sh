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

# Check if required input files exist
if [[ ! -f "$base_dir/MCPS_Freeze_150.GT_hg38.pVCF.rgcpid.QC2.TOPMED_dosages.HDS.pgen" || \
      ! -f "$base_dir/MCPS_Freeze_150.GT_hg38.pVCF.rgcpid.QC2.TOPMED_dosages.HDS.pvar" || \
      ! -f "$output_dir/info8variants.txt" ]]; then
    echo "Error: Required files are missing."
    exit 1
fi

# Create a text file to keep variants with INFO score > 0.8
#awk '$16 > 0.8 {print "chr"$1}' $base_dir/MCPS_Freeze_150.GT_hg38.pVCF.rgcpid.QC2.TOPMED_dosages.snpstats.txt > $chr_out_auto/info8variants.txt


# Phenotype file prefixes and suffixes
phenotype_prefixes=("all-2perc" "all-2.5perc" "all-5perc")
suffixes=("2P-ND" "2.5P-ND" "5P-ND")

# Process phenotype files
for i in "${!phenotype_prefixes[@]}"; do
    prefix="${phenotype_prefixes[$i]}"
    suffix="${suffixes[$i]}"
    
    phenotype_file="$input_dir/phenotypes-ancestry-$suffix.txt"

    # Ensure phenotype files are created correctly
    awk 'BEGIN {print "FID", "IID", "BMI", "AGE", "SEX", "Score_AMR"} {print $2, $2, $3, $4, $5, $25}' \
        "$input_dir/extracted-phenotypes-${prefix}-bmi.txt" > "$phenotype_file"
    
    if [[ ! -f "$phenotype_file" ]]; then
        echo "Error: Phenotype file $phenotype_file was not created."
        exit 1
    fi
    
    echo "Phenotype file $phenotype_file created successfully."

    # Run PLINK2 for each suffix
    $plink2 --pgen $base_dir/MCPS_Freeze_150.GT_hg38.pVCF.rgcpid.QC2.TOPMED_dosages.HDS.pgen \
            --pvar $base_dir/MCPS_Freeze_150.GT_hg38.pVCF.rgcpid.QC2.TOPMED_dosages.HDS.pvar \
            --psam $base_dir/MCPS_Freeze_150.GT_hg38.pVCF.rgcpid.QC2.TOPMED_dosages.HDS.COLLAB.psam \
            --extract $output_dir/info8variants.txt \
            --keep "$phenotype_file" \
            --maf 0.005 \
            --hwe 1e-6 \
            --make-bed \
            --out "$output_dir/topmed-imputed-genotypes-$suffix"

    # Check if PLINK2 ran successfully
    if [[ $? -ne 0 ]]; then
        echo "Error: PLINK2 command failed for $suffix."
        exit 1
    fi

    echo "PLINK2 files for $suffix created successfully."

    # Generate covariates and phenotype files for each suffix
    awk 'BEGIN {print "FID", "IID", "AGE", "SEX", "Score_AMR"} {print $1, $2, $4, $5, $6}' \
        "$phenotype_file" > "$input_dir/covariates_$suffix.cov"

    awk 'BEGIN {print "FID", "IID", "BMI"} {print $1, $2, $3}' \
        "$phenotype_file" > "$input_dir/phenotype_$suffix.pheno"

    echo "Covariates and phenotype files for $suffix created successfully."
done

echo "Job completed successfully."

