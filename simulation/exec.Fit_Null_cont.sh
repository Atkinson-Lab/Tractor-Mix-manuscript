#!/usr/bin/bash
#SBATCH --nodes=1
#SBATCH --partition=mhgcp
#SBATCH --ntasks=1
#SBATCH --job-name=Fit_Null
#SBATCH --output=Fit_Null.%j.out
#SBATCH --mem=32gb
#SBATCH --array=1-6

echo "start"
/opt/R/4.1.1/bin/Rscript Fit_Null_cont.R $SLURM_ARRAY_TASK_ID
echo "finish"
