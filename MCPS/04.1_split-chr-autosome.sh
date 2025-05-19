#!/bin/bash
#SBATCH --job-name=split-bed-navlg
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=05:00:00
#SBATCH -p short
#SBATCH --mail-user=alejandra.vergara@ndph.ox.ac.uk
#SBATCH --array=1-22


# Directories
base_dir=/well/emberson/users/pom143/projects/LAI_TractorMix
data_dir=$base_dir/geno-files/merged-hgdp-1kg-mcps
out_dir=$base_dir/output-files
chr_out_auto=$out_dir/chr_autosomes_sample
plink2=/well/emberson/shared/software/plink2/plink2


chr=${SLURM_ARRAY_TASK_ID}

# Convert to VCF and index
$plink2 --bfile $data_dir/merged_hgdp-1kg-mcps_autosomes --chr $chr --make-bed --out $chr_out_auto/merged_hgdp-1kg-mcps_autosomes_33500_chr$chr

