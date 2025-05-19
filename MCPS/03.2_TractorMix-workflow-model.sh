#!/bin/bash
#SBATCH --job-name=tractormix2-score-model
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=29:00:00
#SBATCH -p short
#SBATCH --mail-user=alejandra.vergara@ndph.ox.ac.uk
#SBATCH --array=21-22

# Record the start time
start_time=$SECONDS

# Activate environment and load necessary modules
module purge
module load Anaconda3/2024.02-1
module load R/4.3.2-gfbf-2023a
module load R-bundle-Bioconductor/3.18-foss-2023a-R-4.3.2


chr=${SLURM_ARRAY_TASK_ID}


Rscript /well/emberson/users/pom143/projects/LAI_TractorMix/src/development/03.2.1_TractorMix-score.R --chr ${chr}
