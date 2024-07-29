#!/usr/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=Power_evaluation_dich
#SBATCH --output=Power_evaluation_dich.%j.out
#SBATCH --mem=32gb
#SBATCH --array=1-99

echo "program start"
/opt/R/4.1.1/bin/Rscript Power_simu_dich.R $SLURM_ARRAY_TASK_ID
echo "program finish"
