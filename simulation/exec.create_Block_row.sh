#!/usr/bin/bash
#SBATCH --nodes=1
#SBATCH --partition=mhgcp
#SBATCH --ntasks=1
#SBATCH --job-name=create_Block_row
#SBATCH --output=create_Block_row.%j.out
#SBATCH --mem=64gb
#SBATCH --array=1-20

echo "start"
/opt/R/4.1.1/bin/Rscript create_Block_row.R $SLURM_ARRAY_TASK_ID
echo "finish"
