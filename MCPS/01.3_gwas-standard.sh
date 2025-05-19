#!/bin/bash
#SBATCH --job-name=topmed-genotype
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --time=20:00:00
#SBATCH -p short
#SBATCH --mail-user=alejandra.vergara@ndph.ox.ac.uk

# Directories
info_dir=/well/emberson/users/pom143/projects/LAI_TractorMix
input_dir=$info_dir/input-files
output_dir=$info_dir/geno-files/mcps
plink2=/well/emberson/shared/software/plink2/plink2

# Phenotype file suffixes
suffixes=("2p" "2.5p" "5p")

# Process phenotype files
for suffix in "${suffixes[@]}"; do
    # Use lowercase suffix for covariate and phenotype files
    cov_file="$input_dir/covariates_$suffix.cov"
    pheno_file="$input_dir/phenotype_$suffix.pheno"
    
    # Convert lowercase suffix to uppercase P for bed, bim, and fam files
    bed_suffix=$(echo "$suffix" | sed 's/p/P-ND/')
    
    bed_file="$output_dir/topmed-imputed-genotypes-$bed_suffix.bed"
    bim_file="$output_dir/topmed-imputed-genotypes-$bed_suffix.bim"
    fam_file="$output_dir/topmed-imputed-genotypes-$bed_suffix.fam"
    
    # Ensure files exist before running PLINK2
    if [[ -f "$cov_file" && -f "$pheno_file" && -f "$bed_file" && -f "$bim_file" && -f "$fam_file" ]]; then
        echo "Processing with suffix: $suffix (genotype suffix: $bed_suffix)"
        
        # Run PLINK2
        $plink2 --bed "$bed_file" \
                --bim "$bim_file" \
                --fam "$fam_file" \
                --covar "$cov_file" \
                --pheno "$pheno_file" \
                --glm \
                --threads 16 \
                --out "$output_dir/plink-$suffix"
                
        if [[ $? -ne 0 ]]; then
            echo "Error: PLINK2 failed for $suffix"
            exit 1
        else
            echo "PLINK2 completed successfully for $suffix"
        fi
    else
        echo "Error: Required files for $suffix not found. Skipping."
    fi
done

echo "All jobs completed."

