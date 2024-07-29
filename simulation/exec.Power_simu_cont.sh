#!/usr/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=Power_evaluation_cont
#SBATCH --output=Power_evaluation_cont.%j.out
#SBATCH --mem=32gb
#SBATCH --array=1-99

echo "program start"
/opt/R/4.1.1/bin/Rscript Power_simu_cont.R $SLURM_ARRAY_TASK_ID
echo "program finish"
