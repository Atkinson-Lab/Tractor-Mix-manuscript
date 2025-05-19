#!/bin/bash
#SBATCH --job-name=grm
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=15G
#SBATCH --time=18:00:00
#SBATCH -p short
#SBATCH --mail-user=alejandra.vergara@ndph.ox.ac.uk
#SBATCH --output=grm-tractormix.out
#SBATCH --error=grm-tractormix.err

# Record the start time
start_time=$SECONDS

# Activate environment and load necessary modules
module purge
module load Anaconda3/2024.02-1
module load R/4.3.2-gfbf-2023a
module load R-bundle-Bioconductor/3.18-foss-2023a-R-4.3.2

# Run the R scripts sequentially with better error handling
set -e  # Exit immediately if a command exits with a non-zero status

# Function to log error and exit
error_exit() {
  echo "Error occurred in script at line $1." >&2
  duration=$(( SECONDS - start_time ))
  echo "The job failed after $((duration / 3600)) hours, $(((duration % 3600) / 60)) minutes, and $((duration % 60)) seconds." >&2
  exit 1
}
trap 'error_exit $LINENO' ERR

# Run the R scripts sequentially
#Rscript /well/emberson/users/pom143/projects/LAI_TractorMix/src/development/03.1.1_TractorMix-pcair_borrar.R
#Rscript /well/emberson/users/pom143/projects/LAI_TractorMix/src/development/03.1.2_TractorMix-pcrel-vf.R
Rscript /well/emberson/users/pom143/projects/LAI_TractorMix/src/development/03.1.3_TractorMix-grm-model.R

# Calculate and print the duration for the job
duration=$(( SECONDS - start_time ))
echo "The job took $((duration / 3600)) hours, $(((duration % 3600) / 60)) minutes, and $((duration % 60)) seconds."

