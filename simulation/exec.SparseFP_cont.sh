#!/usr/bin/bash
#SBATCH --nodes=1
#SBATCH --partition=mhgcp
#SBATCH --ntasks=1
#SBATCH --job-name=SparseFP_cont
#SBATCH --output=SparseFP_cont.%j.out
#SBATCH --mem=32gb
#SBATCH --array=1-20

echo "start"
/opt/R/4.1.1/bin/Rscript SparseFP_cont.R $SLURM_ARRAY_TASK_ID
echo "finish"
