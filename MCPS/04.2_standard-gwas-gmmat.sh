#!/bin/bash
#SBATCH --job-name=standard-gwas-navlg
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=12G
#SBATCH --time=29:00:00
#SBATCH -p short
#SBATCH --mail-user=alejandra.vergara@ndph.ox.ac.uk
#SBATCH --array=1-22

# Activate environment and load necessary modules
module purge
module load Anaconda3/2024.02-1
module load R/4.3.2-gfbf-2023a
module load R-bundle-Bioconductor/3.18-foss-2023a-R-4.3.2

# Directories
tractor_dir=/well/emberson/users/pom143/projects/LAI_TractorMix

chr=${SLURM_ARRAY_TASK_ID}

# Run the first R script and, if successful, run the second R script
Rscript $tractor_dir/src/development/04.2.1_standard-gwas-full-parallel.R --chr ${chr}
#Rscript $tractor_dir/src/development/04.2.2_standard_gwas_sparse.R --chr ${chr}
