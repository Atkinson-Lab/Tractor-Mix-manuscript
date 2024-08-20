#!/usr/bin/bash
#SBATCH --nodes=1
#SBATCH --partition=mhgcp
#SBATCH --nodelist=mhgcp-c00,mhgcp-c01,mhgcp-c02,mhgcp-m00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=100:00:00
#SBATCH --job-name=Tractor-Mix
#SBATCH --output=Tractor-Mix.%j.out
#SBATCH --mem=32gb
#SBATCH --array=1-22

echo "start TractorMix..."
/opt/R/4.1.1/bin/Rscript TC_LDL_SCA_TractorMix.R $SLURM_ARRAY_TASK_ID
echo "finish TractorMix..."
